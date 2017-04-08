// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sys/types.h>
#include <sys/fcntl.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <unistd.h>

#if !defined(__STDC_NO_THREADS__) || __STDC_NO_THREADS__
#define USE_PTHREADS 1
#else
#define USE_PTHREADS 0
#endif

#if USE_PTHREADS
#include <pthread.h>
#else
#include <threads.h>
#endif

#include "mpi.h"
#include "sodium.h"
#include "core/polymec.h"
#include "core/options.h"
#include "core/array.h"
#include "core/unordered_map.h"
#include "core/engine_server.h"

// We only do this stuff on MPI rank 0.
static int _rank = -1;

// This creates the resource directory and populates it with needed items.
static char _resource_dir[FILENAME_MAX+1];
static void set_up_resource_dir(const char* resource_dir)
{
  ASSERT(resource_dir != NULL);
  create_directory(resource_dir, S_IRUSR | S_IRGRP | S_IWUSR | S_IXUSR);
  strncpy(_resource_dir, resource_dir, FILENAME_MAX);
}

//------------------------------------------------------------------------
//                             Key Database
//------------------------------------------------------------------------

static const char* _keys_dirname = ".keys";
static const char* _s_pub_key_filename = "engine.pub";
static const char* _s_priv_key_filename = "engine";
static const char* _c_pub_key_filename = "clients";
static const char* _keys_lock_filename = ".keys.lock";

typedef struct
{
  // Server key pair.
  size_t s_priv_key_len, s_pub_key_len;
  uint8_t *s_priv_key, *s_pub_key;

  // Client public keys.
  string_ptr_unordered_map_t* c_keys;
  bool updated;

  // Lock file descriptor.
  int lock;
} keydb_t;

typedef struct
{
  char* username;
  char* password;
  uint8_t* key;
  size_t key_len;
} keydb_entry_t;

static keydb_entry_t* keydb_entry_new(const char* username,
                                      const char* password,
                                      size_t key_len,
                                      uint8_t key[key_len])
{
  keydb_entry_t* entry = polymec_malloc(sizeof(keydb_entry_t));
  entry->username = string_dup(username);
  entry->password = string_dup(password);
  entry->key_len = key_len;
  entry->key = polymec_malloc(sizeof(uint8_t) * key_len);
  memcpy(entry->key, key, sizeof(uint8_t) * key_len);
  return entry;
}

static void keydb_entry_free(void* e)
{
  keydb_entry_t* entry = e;
  string_free(entry->username);
  string_free(entry->password);
  polymec_free(entry->key);
  polymec_free(entry);
}

static keydb_t* keydb_open()
{
  ASSERT(_rank == 0);
  ASSERT(directory_exists(_resource_dir));

  // Create a lock file and lock it.
  char lock_file[FILENAME_MAX+1];
  snprintf(lock_file, FILENAME_MAX, "%s/%s", _resource_dir, _keys_lock_filename);
  int lock = open(lock_file, O_WRONLY);
  struct flock fl;
  fl.l_type = F_WRLCK;
  fl.l_whence = SEEK_SET;
  fl.l_start = 0;
  fl.l_len = 0;
  fl.l_pid = getpid();
  fcntl(lock, F_SETLKW, &fl);

  // Make sure we have a directory to store our keys.
  char keys_dir[FILENAME_MAX+1];
  snprintf(keys_dir, FILENAME_MAX, "%s/%s", _resource_dir, _keys_dirname);
  if (!directory_exists(keys_dir))
    create_directory(keys_dir, S_IRUSR | S_IRGRP | S_IWUSR | S_IXUSR);

  // Do we already have a key pair for this server? If not, create one.
  char s_pub_key_file[FILENAME_MAX+1], s_priv_key_file[FILENAME_MAX+1];
  snprintf(s_pub_key_file, FILENAME_MAX, "%s/%s", keys_dir, _s_pub_key_filename);
  snprintf(s_priv_key_file, FILENAME_MAX, "%s/%s", keys_dir, _s_priv_key_filename);
  if (!file_exists(s_pub_key_file) || !file_exists(s_priv_key_file))
  {
  }

  // Now check to see whether we have a client key database.
  char c_pub_key_file[FILENAME_MAX+1];
  snprintf(c_pub_key_file, FILENAME_MAX, "%s/%s", keys_dir, _c_pub_key_filename);
  if (!file_exists(c_pub_key_file))
  {
  }

  // Load the keys into memory.
  keydb_t* keys = polymec_malloc(sizeof(keydb_t));
  keys->c_keys = string_ptr_unordered_map_new();
  keys->updated = false;
  keys->lock = lock;
  return keys;
}

static void keydb_close(keydb_t* keys)
{
  // Write out our client key database if we have updates.
  if (keys->updated)
  {
  }

  // Release the lock file.
  struct flock fl;
  fl.l_type = F_UNLCK;
  fl.l_whence = SEEK_SET;
  fl.l_start = 0;
  fl.l_len = 0;
  fl.l_pid = getpid();
  fcntl(keys->lock, F_SETLKW, &fl);
  close(keys->lock);

  // Delete it if we can.
  char lock_file[FILENAME_MAX+1];
  snprintf(lock_file, FILENAME_MAX, "%s/%s", _resource_dir, _keys_lock_filename);
  if (file_exists(lock_file))
    remove(lock_file);

  polymec_free(keys->s_priv_key);
  polymec_free(keys->s_pub_key);
  string_ptr_unordered_map_free(keys->c_keys);
  polymec_free(keys);
}

static uint8_t* keydb_server_pub_key(keydb_t* keys, size_t* key_len)
{
  return keys->s_pub_key;
}

static uint8_t* keydb_server_priv_key(keydb_t* keys, size_t* key_len)
{
  return keys->s_priv_key;
}

static uint8_t* keydb_client_key(keydb_t* keys, const char* username, size_t* key_len)
{
  keydb_entry_t** e_ptr = (keydb_entry_t**)string_ptr_unordered_map_get(keys->c_keys, (char*)username);
  if (e_ptr != NULL) 
  {
    *key_len = (*e_ptr)->key_len;
    return (*e_ptr)->key;
  }
  else 
    return NULL;
}

static void keydb_insert_user(keydb_t* keys,
                              const char* username,
                              const char* password,
                              size_t key_len,
                              uint8_t public_key[key_len])
{
  keydb_entry_t* e = keydb_entry_new(username, password, key_len, public_key);
  string_ptr_unordered_map_insert_with_kv_dtors(keys->c_keys, (char*)username, e, 
                                                string_free,
                                                keydb_entry_free);
  keys->updated = true;
}

static void keydb_delete_user(keydb_t* keys,
                              const char* username)
{
  string_ptr_unordered_map_delete(keys->c_keys, (char*)username);
  keys->updated = true;
}

static bool keydb_contains_user(keydb_t* keys,
                                const char* username)
{
  return string_ptr_unordered_map_contains(keys->c_keys, (char*)username);
}

//------------------------------------------------------------------------
//                            Server machinery
//------------------------------------------------------------------------

// Server address and socket.
struct sockaddr_in _server_addr;
static int _server = -1;

// Server identity file name.
char _id_file[FILENAME_MAX+1];

// Queue of outgoing messages from server -> client.
ptr_array_t* _outgoing_mesgs = NULL;

// Listener thread.
#if USE_PTHREADS
pthread_attr_t _listener_attr;
pthread_t _listener;
#else
thrd_t _listener;
#endif
bool _shutting_down = false;

// This authenticates a client that has connected, returning true if successful
// and false otherwise.
static bool authenticate_client(int client)
{
  ASSERT(_rank == 0);

  // Open our key database.
  keydb_t* keys = keydb_open();

  // Get the username and password from the client.
  char* username;
  char* password;
  if (!keydb_contains_user(keys, (const char*)username))
    return false;

  // Read our private key.
  size_t s_priv_key_len;
  uint8_t *s_priv_key = keydb_server_priv_key(keys, &s_priv_key_len);

  // Retrieve the client's public key from our database.
  size_t c_pub_key_len;
  uint8_t* c_pub_key = keydb_client_key(keys, username, &c_pub_key_len);

  // Close the key database.
  keydb_close(keys);

  return true;
}

// This handles client requests until the client disconnects.
static void handle_client_requests(int client)
{
  while (true)
  {
  }
}

// This accepts connections from clients.
#if USE_PTHREADS
static void* accept_connections(void* context)
#else
static int accept_connections(void* context)
#endif
{
  ASSERT(_server != -1);
  _shutting_down = false;

  // This thread should be cancelable.
#if USE_PTHREADS
  int old;
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &old);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, &old);
#else
#endif

  // Listen on the server socket.
  listen(_server, 5);

  while (true)
  {
    // Are we shutting down?
    if (_shutting_down) 
      break;

    // Accept a new connection.
    struct sockaddr_in client_addr;
    socklen_t client_len = sizeof(client_addr);
    int client = accept(_server, (struct sockaddr*)&client_addr, &client_len);
    if (client < 0)
      log_info("accept_connections: Could not accept client connection.");

    // Authenticate.
    if (authenticate_client(client))
    {
      // Handle the client's requests.
      handle_client_requests(client);
    }

    // Close the client socket.
    close(client);
  }

  // Get out.
#if USE_PTHREADS
  pthread_exit(NULL);
  return NULL;
#else
  thrd_exit(0);
  return 0;
#endif
}

static void stop_engine_server()
{
  if (engine_server_is_running())
    engine_server_stop();
}

//------------------------------------------------------------------------
//                              Server API
//------------------------------------------------------------------------

bool engine_server_is_running()
{
  return (_server != -1);
}

void engine_server_start(const char* resource_dir)
{
  ASSERT(_server == -1);
  ASSERT(resource_dir != NULL);

  // Find our MPI rank.
  if (_rank == -1)
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

  if (_rank == 0)
  {
    if (!directory_exists(resource_dir))
      set_up_resource_dir(resource_dir);

    options_t* opts = options_argv();

    // Engine server port number.
    char* port_str = options_value(opts, "engine_port");

    if ((port_str != NULL) && 
        string_is_number((const char*)port_str))
    {
      // Get the port number.
      int port = atoi((const char*)port_str);

      // Try to listen on that port.
      _server = socket(AF_INET, SOCK_STREAM, 0);
      if (_server != -1)
      {
        // Set up our server's address.
        memset(&_server_addr, 0, sizeof(struct sockaddr_in));
        _server_addr.sin_family = AF_INET;
        _server_addr.sin_addr.s_addr = INADDR_ANY;
        _server_addr.sin_port = htons(port);

        // Try to bind the address so we can listen on it.
        int status = bind(_server, (struct sockaddr*)&_server_addr, sizeof(_server_addr));
        if (status >= 0)
        {
          // Fire up sodium and register our shutdown function.
          if (sodium_init() != -1)
          {
            polymec_atexit(stop_engine_server);

            // Spawn our listener thread.
#if USE_PTHREADS
            pthread_attr_init(&_listener_attr);
            pthread_create(&_listener, &_listener_attr, accept_connections, NULL);
#else
            thrd_create(&_listener, accept_connections, NULL);
#endif
          }
          else
            log_info("engine_server_start: Couldn't initialize encryption library. Not starting server.");
        }
        else
          log_info("engine_server_start: Couldn't bind port %d for the compute engine server. Not starting server.", port);
      }
      else
        log_info("engine_server_start: Couldn't open a socket for the compute engine server. Not starting server.");
    }
  }

  // Broadcast the value of _server to other ranks.
  MPI_Bcast(&_server, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void engine_server_stop()
{
  ASSERT(_rank != -1);
  ASSERT(_server != -1);

  if (_rank == 0)
  {
    _shutting_down = true;

    // Wait for the thread to join.
#if USE_PTHREADS
    void* retval;
    pthread_join(_listener, &retval);
#else
    int retval;
    thrd_join(_listener, &retval);
#endif
  }
}

void engine_server_halt()
{
  ASSERT(_rank != -1);
  ASSERT(_server != -1);
  if (_rank == 0)
  {
    shutdown(_server, SHUT_RDWR);
#if USE_PTHREADS
    pthread_cancel(_listener);
#else
    engine_server_stop(); // C11 doesn't have noncooperative cancellation. :-/
#endif
  }
}

void engine_server_insert_user(const char* username, 
                               const char* password,
                               size_t key_len,
                               uint8_t public_key[key_len])
{
  ASSERT(_server != -1);
  keydb_t* keys = keydb_open();
  keydb_insert_user(keys, username, password, key_len, public_key);
  keydb_close(keys);
}

void engine_server_delete_user(const char* username)
{
  ASSERT(_server != -1);
  keydb_t* keys = keydb_open();
  keydb_delete_user(keys, username);
  keydb_close(keys);
}

uint8_t* engine_server_public_key(size_t* key_len)
{
  ASSERT(_server != -1);
  keydb_t* keys = keydb_open();
  uint8_t* key = keydb_server_pub_key(keys, key_len);
  keydb_close(keys);
  return key;
}

