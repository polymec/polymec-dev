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
#include "core/slist.h"
#include "core/array.h"
#include "core/serializer.h"
#include "core/unordered_map.h"
#include "core/engine_server.h"

// We only do this stuff on MPI rank 0.
static int _rank = -1;

//------------------------------------------------------------------------
//                             Key Database
//------------------------------------------------------------------------

// This creates the resource directory and populates it with needed items.
static char _resource_dir[FILENAME_MAX+1];
static void set_up_resource_dir(const char* resource_dir)
{
  ASSERT(resource_dir != NULL);
  create_directory(resource_dir, S_IRUSR | S_IRGRP | S_IWUSR | S_IXUSR);
  strncpy(_resource_dir, resource_dir, FILENAME_MAX);
}

static const char* _keys_dirname = ".keys";
static const char* _s_pub_key_filename = "engine.pub";
static const char* _s_priv_key_filename = "engine";
static const char* _c_pub_key_filename = "clients";
static const char* _keys_lock_filename = ".keys.lock";

static char _keys_dir[FILENAME_MAX+1];

typedef struct
{
  // Server key pair.
  uint8_t s_priv_key[crypto_box_SECRETKEYBYTES], s_pub_key[crypto_box_PUBLICKEYBYTES];

  // Client public keys.
  string_ptr_unordered_map_t* c_keys;
  bool updated;

  // Lock file descriptor.
  int lock;
} keydb_t;

typedef struct
{
  char* username;
  char* hashed_password;
  uint8_t* key;
} keydb_entry_t;

static keydb_entry_t* keydb_entry_new(const char* username,
                                      const char* hashed_password,
                                      uint8_t* key)
{
  keydb_entry_t* entry = polymec_malloc(sizeof(keydb_entry_t));
  entry->username = string_dup(username);
  entry->hashed_password = string_dup(hashed_password);
  entry->key = polymec_malloc(sizeof(uint8_t) * crypto_box_PUBLICKEYBYTES);
  memcpy(entry->key, key, sizeof(uint8_t) * crypto_box_PUBLICKEYBYTES);
  return entry;
}

static void keydb_entry_free(void* e)
{
  keydb_entry_t* entry = e;
  string_free(entry->username);
  string_free(entry->hashed_password);
  polymec_free(entry->key);
  polymec_free(entry);
}

static bool keydb_insert_user(keydb_t* keys,
                              const char* username,
                              const char* password,
                              uint8_t* public_key)
{
  bool inserted = false;
  if (!string_ptr_unordered_map_contains(keys->c_keys, (char*)username))
  {
    keydb_entry_t* e = keydb_entry_new(username, password, public_key);
    string_ptr_unordered_map_insert_with_kv_dtors(keys->c_keys, (char*)username, e, 
                                                  string_free,
                                                  keydb_entry_free);
    keys->updated = true;
    inserted = true;
  }
  return inserted;
}

static void keydb_gen_server_keys(keydb_t* keys)
{
  // Open the server key files for writing.
  char s_pub_key_file[FILENAME_MAX+1], s_priv_key_file[FILENAME_MAX+1];
  snprintf(s_pub_key_file, FILENAME_MAX, "%s/%s", _keys_dir, _s_pub_key_filename);
  snprintf(s_priv_key_file, FILENAME_MAX, "%s/%s", _keys_dir, _s_priv_key_filename);
  if (!file_exists(s_pub_key_file) || !file_exists(s_priv_key_file))
  {
    FILE* s_pub = fopen(s_pub_key_file, "w");
    fclose(s_pub);

    FILE* s_priv = fopen(s_priv_key_file, "w");
    fclose(s_priv);
  }
}

static void keydb_write_client_keys(keydb_t* keys, const char* key_file)
{
  // Open the client key database for writing.
  char c_pub_key_file[FILENAME_MAX+1];
  snprintf(c_pub_key_file, FILENAME_MAX, "%s/%s", _keys_dir, key_file);
  FILE* c_keys = fopen(c_pub_key_file, "w");

  // Loop through the client entries and write them out.
  int pos = 0;
  char* username;
  keydb_entry_t* entry;
  while (string_ptr_unordered_map_next(keys->c_keys, &pos, &username, (void**)(&entry)))
  {
    // Write the username and a : delimiter.
    fprintf(c_keys, "%s:", entry->username);

    // Write the hashed password and a : delimiter.
    fprintf(c_keys, "%s:", entry->hashed_password);

    // Write the public key.
    for (size_t i = 0; i < crypto_box_PUBLICKEYBYTES; ++i)
      fprintf(c_keys, "%02u", entry->key[i]);

    // Newline.
    fprintf(c_keys, "\n");
  }

  // Close the file.
  fclose(c_keys);
}

static void keydb_read_client_keys(keydb_t* keys, const char* key_file)
{
  // Open the client key database for writing.
  char c_pub_key_file[FILENAME_MAX+1];
  snprintf(c_pub_key_file, FILENAME_MAX, "%s/%s", _keys_dir, key_file);
  FILE* c_keys = fopen(c_pub_key_file, "r");

  // Loop through the client entries and read them.
  int n = 0;
  while (n != EOF)
  {
    // Read the username and a : delimiter.
    char username[32+1];
    n = fscanf(c_keys, "%32s:", username);
    if (n < 0)
    {
      log_info("keydb_read_client_keys: Invalid client username read.");
      continue;
    }
    username[n] = '\0';

    // Read the hashed password and a : delimiter.
    char hashed_pw[32+1];
    n = fscanf(c_keys, "%32s:", hashed_pw);
    if (n < 0)
    {
      log_info("keydb_read_client_keys: Invalid client password read.");
      continue;
    }
    hashed_pw[n] = '\0';

    // Read the public key.
    uint8_t key[crypto_box_PUBLICKEYBYTES];
    for (size_t i = 0; i < crypto_box_PUBLICKEYBYTES; ++i)
    {
      n = fscanf(c_keys, "%02u", &key[i]);
      if (n != 2)
      {
        log_info("keydb_read_client_keys: Invalid client key read.");
        continue;
      }
    }

    // Read the newline.
    n = fscanf(c_keys, "\n");
    if (n != 1)
    {
      log_info("keydb_read_client_keys: Invalid line ending read.");
      continue;
    }

    // Insert the entry.
    keydb_insert_user(keys, username, hashed_pw, key);
  }

  // Close the file.
  fclose(c_keys);
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
  snprintf(_keys_dir, FILENAME_MAX, "%s/%s", _resource_dir, _keys_dirname);
  if (!directory_exists(_keys_dir))
    create_directory(_keys_dir, S_IRUSR | S_IRGRP | S_IWUSR | S_IXUSR);

  // Create the database. 
  keydb_t* keys = polymec_malloc(sizeof(keydb_t));
  keys->c_keys = string_ptr_unordered_map_new();
  keys->updated = false;
  keys->lock = lock;

  // Do we already have a key pair for this server? If not, create one.
  keydb_gen_server_keys(keys);

  // Now check to see whether we have a client key database.
  char c_pub_key_file[FILENAME_MAX+1];
  snprintf(c_pub_key_file, FILENAME_MAX, "%s/%s", _keys_dir, _c_pub_key_filename);
  if (!file_exists(c_pub_key_file))
    keydb_write_client_keys(keys, c_pub_key_file);
  else
    keydb_read_client_keys(keys, c_pub_key_file);

  return keys;
}

static void keydb_close(keydb_t* keys)
{
  // Write out our client key database if we have updates.
  if (keys->updated)
  {
    char c_pub_key_file[FILENAME_MAX+1];
    snprintf(c_pub_key_file, FILENAME_MAX, "%s/%s", _keys_dir, _c_pub_key_filename);
    keydb_write_client_keys(keys, c_pub_key_file);
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

static uint8_t* keydb_server_pub_key(keydb_t* keys)
{
  return keys->s_pub_key;
}

//static uint8_t* keydb_server_priv_key(keydb_t* keys)
//{
//  return keys->s_priv_key;
//}

static uint8_t* keydb_client_key(keydb_t* keys, 
                                 const char* username,
                                 const char* hashed_password)
{
  keydb_entry_t** e_ptr = (keydb_entry_t**)string_ptr_unordered_map_get(keys->c_keys, (char*)username);
  if (e_ptr != NULL) 
  {
    if (strcmp((*e_ptr)->hashed_password, hashed_password) == 0)
      return (*e_ptr)->key;
    else
      return NULL;
  }
  else 
    return NULL;
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
ptr_slist_t* _queue = NULL;
#if USE_PTHREADS
pthread_mutex_t _queue_lock;
pthread_cond_t _queue_cond;
#else
mtx_t _queue_lock;
cnd_t _queue_cond;
#endif

// Listener thread.
#if USE_PTHREADS
pthread_attr_t _listener_attr;
pthread_t _listener;
#else
thrd_t _listener;
#endif
bool _shutting_down = false;

static uint8_t* authenticate_client(int client)
{
  // Get the username and the hashed password for our client.

  // Request username.
  write(_server, "OHAI", sizeof(char)*4);

  // Receive username.
  char username[32+1];
  ssize_t n = read(client, username, 32);
  if (n < 0) 
  {
    log_info("authenticate_client: Could not read username.");
    return NULL;
  }
  username[n] = '\0';

  // Request hashed password.
  write(_server, "A/S/L", sizeof(char)*5);
  char hashed_pw[32+1];
  n = read(client, hashed_pw, 32);
  if (n < 0) 
  {
    log_info("authenticate_client: Could not read password.");
    return NULL;
  }
  hashed_pw[n] = '\0';

  // Get the key.
  uint8_t* key = NULL;
  {
    keydb_t* keys = keydb_open();
    uint8_t* k = keydb_client_key(keys, username, hashed_pw);
    if (k != NULL)
    {
      key = polymec_malloc(sizeof(uint8_t) * crypto_box_PUBLICKEYBYTES);
      memcpy(key, k, sizeof(uint8_t) * crypto_box_PUBLICKEYBYTES);
    }
    else
      log_info("authenticate_client: Authentication failed.");
    keydb_close(keys);
  }
  return key;
}

// Methods for managing access to the message queue.
static inline void lock_queue()
{
#if USE_PTHREADS
  pthread_mutex_lock(&_queue_lock);
#else
  mtx_lock(&_queue_lock);
#endif
}

static inline void unlock_queue()
{
#if USE_PTHREADS
  pthread_mutex_unlock(&_queue_lock);
#else
  mtx_unlock(&_queue_lock);
#endif
}

static inline void wait_for_messages()
{
  while (!_shutting_down)
#if USE_PTHREADS
    pthread_cond_wait(&_queue_cond, &_queue_lock);
#else
    cnd_wait(&_queue_cond, &_queue_lock);
#endif
}

// This handles client requests until the client disconnects.
static void handle_client_requests(int client, uint8_t* client_key)
{
  while (true)
  {
    // Lock the message queue and wait for messages.
    lock_queue();
    wait_for_messages();

    // If we're asked to shut down, do so quietly.
    if (_shutting_down)
    {
      unlock_queue();
      break;
    }

    if (!ptr_slist_empty(_queue))
    {
      // Grab the next message.
      void (*destroy_msg)(void*);
      byte_array_t* msg = ptr_slist_pop(_queue, &destroy_msg);

      // Encrypt the message using our private key and the client's
      // public key.

      // Send the encrypted message.

      // Destroy the message on our end.
      destroy_msg(msg);
    }

    unlock_queue();
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
    uint8_t* client_key = authenticate_client(client);
    if (client_key != NULL)
    {
      // Handle the client's requests.
      handle_client_requests(client, client_key);

      // Dispose of the key.
      polymec_free(client_key);
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
            polymec_atexit(engine_server_stop);

            // Set stuff up.
#if USE_PTHREADS
            pthread_mutex_init(&_queue_lock, NULL);
            pthread_cond_init(&_queue_cond, NULL);
#else
            mtx_init(&_queue_lock, mtx_plain);
            cnd_init(&_queue_cond);
#endif

            // Lock the queue till we're ready to go.
            lock_queue();

            // Spawn our listener thread.
#if USE_PTHREADS
            pthread_attr_init(&_listener_attr);
            pthread_create(&_listener, &_listener_attr, accept_connections, NULL);
#else
            thrd_create(&_listener, accept_connections, NULL);
#endif
            // Okay, do it.
            unlock_queue();
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

static void cleanup()
{
  // Kill the queue.
  if (_queue != NULL)
  {
    ptr_slist_free(_queue);
    _queue = NULL;

#if USE_PTHREADS
    pthread_mutex_destroy(&_queue_lock);
    pthread_cond_destroy(&_queue_cond);
    pthread_attr_destroy(&_listener_attr);
#else
    mtx_destroy(&_queue_lock);
    cnd_destroy(&_queue_cond);
#endif
  }
  _rank = -1;
}

void engine_server_stop()
{
  if ((_rank == 0) && (_server != -1))
  {
    // Hey, listener! We're shutting down.
    lock_queue();
    _shutting_down = true;
#if USE_PTHREADS
    pthread_cond_signal(&_queue_cond);
#else
    cnd_signal(&_queue_cond);
#endif
    unlock_queue();

    // Wait for the thread to join.
#if USE_PTHREADS
    void* retval;
    pthread_join(_listener, &retval);
#else
    int retval;
    thrd_join(_listener, &retval);
#endif
  }

  cleanup();
}

void engine_server_halt()
{
  if ((_rank == 0) && (_server != -1))
  {
    shutdown(_server, SHUT_RDWR);
#if USE_PTHREADS
    pthread_cancel(_listener);
#else
    engine_server_stop(); // C11 doesn't have noncooperative cancellation. :-/
#endif
  }
  cleanup();
}

char* engine_server_hash_password(const char* plain_text_password)
{
  char* hashed_pw = polymec_malloc(sizeof(char) * crypto_pwhash_STRBYTES);
  int stat = crypto_pwhash_str(hashed_pw, plain_text_password, strlen(plain_text_password),
                               crypto_pwhash_OPSLIMIT_SENSITIVE, 
                               crypto_pwhash_MEMLIMIT_SENSITIVE);
  ASSERT(stat != 0); // Out of memory?
  ASSERT(crypto_pwhash_str_verify(hashed_pw, plain_text_password, strlen(plain_text_password)) == 0);

  return hashed_pw;
}

bool engine_server_insert_user(const char* username, 
                               const char* password,
                               uint8_t* public_key)
{
  bool inserted = false;
  if (engine_server_is_running())
  {
    keydb_t* keys = keydb_open();
    if (!keydb_contains_user(keys, username))
      inserted = keydb_insert_user(keys, username, password, public_key);
    keydb_close(keys);
  }
  return inserted;
}

void engine_server_delete_user(const char* username)
{
  if (engine_server_is_running())
  {
    keydb_t* keys = keydb_open();
    keydb_delete_user(keys, username);
    keydb_close(keys);
  }
}

uint8_t* engine_server_public_key()
{
  uint8_t* key = NULL;
  if (engine_server_is_running())
  {
    keydb_t* keys = keydb_open();
    uint8_t* k = keydb_server_pub_key(keys);
    if (k != NULL)
    {
      key = polymec_malloc(sizeof(uint8_t) * crypto_box_PUBLICKEYBYTES);
      memcpy(key, k, sizeof(uint8_t) * crypto_box_PUBLICKEYBYTES);
    }
    keydb_close(keys);
  }
  return key;
}

void engine_server_push_command(const char* command)
{
  ASSERT(sizeof(char) == sizeof(uint8_t));
  if ((_rank == 0) && engine_server_is_running())
  {
    size_t len = strlen(command);
    ASSERT(len <= 128);
    byte_array_t* data = byte_array_new();
    byte_array_resize(data, len);
    memcpy(data->data, command, len*sizeof(char));
    engine_server_push_data(data);
  }
}

void engine_server_push_data(byte_array_t* data)
{
  if ((_rank == 0) && engine_server_is_running())
  {
    // Acquire the lock for the message queue.
    lock_queue();

    if (_queue == NULL)
      _queue = ptr_slist_new();

    ptr_slist_append_with_dtor(_queue, data, DTOR(byte_array_free));

    // Let the listener thread know there's a new message.
#if USE_PTHREADS
    pthread_cond_signal(&_queue_cond);
#else
    cnd_signal(&_queue_cond);
#endif

    // Unlock the queue.
    unlock_queue();
  }
}

void engine_server_push_real_array(int rank, size_t* shape, real_t* array)
{
  byte_array_t* bytes = byte_array_new();
  size_t offset = 0;
  byte_array_write_ints(bytes, 1, &rank, &offset);
  byte_array_write_size_ts(bytes, (size_t)rank, shape, &offset);
  size_t n = 1;
  for (int i = 0; i < rank; ++i)
    n *= shape[i];
  byte_array_write_real_ts(bytes, n, array, &offset);
  polymec_free(array);
}

void engine_server_push_complex_array(int rank, size_t* shape, complex_t* array)
{
  byte_array_t* bytes = byte_array_new();
  size_t offset = 0;
  byte_array_write_ints(bytes, 1, &rank, &offset);
  byte_array_write_size_ts(bytes, (size_t)rank, shape, &offset);
  size_t n = 1;
  for (int i = 0; i < rank; ++i)
    n *= shape[i];
  byte_array_write_complex_ts(bytes, n, array, &offset);
  polymec_free(array);
}
