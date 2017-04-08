// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>

#if !defined(__STDC_NO_THREADS__) || __STDC_NO_THREADS__
#include <pthread.h>
#else
#include <threads.h>
#endif

#include "sodium.h"
#include "core/polymec.h"
#include "core/options.h"
#include "core/array.h"

// This creates the resource directory and populates it with needed items.
char _resource_dir[FILENAME_MAX+1];
static void set_up_resource_dir(const char* resource_dir)
{
  if (resource_dir != NULL)
  {
    if (!directory_exists(resource_dir))
    {
      bool created = create_directory(resource_dir, 
          S_IRUSR | S_IRGRP | S_IWUSR | S_IXUSR);
      if (!created)
        polymec_error("Cannot create resource directory %s.", resource_dir);
    }
    strncpy(_resource_dir, resource_dir, FILENAME_MAX);
  }
  else
    _resource_dir[0] = '\0';
}

// Server address and socket.
struct sockaddr_in _server_addr;
static int _server = -1;

// Server identity file name.
char _id_file[FILENAME_MAX+1];

// Queue of outgoing messages from server -> client.
ptr_array_t* _outgoing_mesgs = NULL;

// Listener thread.
#if !defined(__STDC_NO_THREADS__) || __STDC_NO_THREADS__
pthread_attr_t _listener_attr;
pthread_t _listener;
#else
thrd_t _listener;
#endif
bool _shutting_down = false;

// This shuts down the compute engine server.
static void server_stop(void)
{
  _shutting_down = true;
}

// This authenticates a client that has connected, returning 0 if successful
// and -1 otherwise.
static int authenticate_client(int client)
{
  // Read our public and private keys from our key file.
  size_t server_pub_key_len, server_priv_key_len;
  uint8_t server_pub_key[server_pub_key_len], server_priv_key[server_priv_key_len];

  // Retrieve the client's public key from our database.
  size_t client_pub_key_len;
  uint8_t client_pub_key[client_pub_key_len];
}

// This handles client requests until the client disconnects.
static void handle_client_requests(int client)
{
  while (true)
  {
  }
}

// This accepts connections from clients.
#if __STDC_NO_THREADS__
static void* accept_connections(void* context)
#else
static int accept_connections(void* context)
#endif
{
  ASSERT(_server != -1);
  _shutting_down = false;

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
    if (authenticate_client(client) == 0)
    {
      // Handle the client's requests.
      handle_client_requests(client);
    }

    // Close the client socket.
    close(client);
  }

  // Get out.
#if !defined(__STDC_NO_THREADS__) || __STDC_NO_THREADS__
  pthread_exit(NULL);
  return NULL;
#else
  thrd_exit(0);
  return 0;
#endif
}

// This generates public and private keys for the server if requested.
const char* _key_dirname = ".keys";
const char* _pub_key_filename = "engine.pub";
const char* _priv_key_filename = "engine";
void polymec_server_keygen(void);
void polymec_server_keygen()
{
  const char* resource_dir = polymec_resource_dir();
  if (resource_dir == NULL)
    polymec_error("polymec_server_keygen: Can't generate keys for this executeable (no resource directory).");
  ASSERT(directory_exists(resource_dir));

  char key_dir[FILENAME_MAX+1];
  snprintf(key_dir, FILENAME_MAX, "%s/%s", resource_dir, _key_dirname);
}

// This starts the compute engine server if requested.
void polymec_server_start(void)
{
  ASSERT(_server == -1);

  options_t* opts = options_argv();
  
  // Compute server port number.
  char* port_str = options_value(opts, "server_port");

  // If we can't bind the port, do we fail?
  char* fail_str = options_value(opts, "server_required"); 
  bool fail = (string_is_boolean(fail_str) && string_as_boolean(fail_str));

  // Do we have a key pair for this server? 
  bool has_keys = false;
  const char* resource_dir = polymec_resource_dir();
  if (resource_dir != NULL)
  {
    char key_dir[FILENAME_MAX+1];
    snprintf(key_dir, FILENAME_MAX, "%s/%s", resource_dir, _key_dirname);
    if (directory_exists(key_dir))
    {
      char pub_key[FILENAME_MAX+1], priv_key[FILENAME_MAX+1];
      snprintf(pub_key, FILENAME_MAX, "%s/%s", key_dir, _pub_key_filename);
      snprintf(priv_key, FILENAME_MAX, "%s/%s", key_dir, _priv_key_filename);
      if (file_exists(pub_key) && file_exists(priv_key))
        has_keys = true;
    }
  }
  if (!has_keys)
  {
    if (fail)
      polymec_error("polymec_server_start: No key pair found for encryption. Please run with the --keygen option");
    else
      log_info("polymec_server_start: No key pair found for encryption. Won't start server.");
  }
  else
  {
    if (!file_exists(id_file))
    {
      if (fail)
        polymec_error("polymec_server_start: Can't find encryption identity file: %s", id_file);
      else
        log_info("polymec_server_start: Can't find encryption identity file: %s. Won't start server.", id_file);
    }

    // Copy the identity file into place.
    strncpy(_id_file, id_file, FILENAME_MAX);
  }

  if ((port_str != NULL) && 
      string_is_number((const char*)port_str) && 
      (id_file != NULL))
  {
    // Get the port number.
    int port = atoi((const char*)port_str);

    // Try to listen on that port.
    _server = socket(AF_INET, SOCK_STREAM, 0);
    if (_server != -1)
    {
      // Set up our server's address.
      memset(&_server_addr, 0, sizeof(struct sockaddr_in))
      _server_addr.sin_family = AF_INET;
      _server_addr.sin_addr.s_addr = INADDR_ANY;
      _server_addr.sin_port = htons(port);

      // Try to bind the address so we can listen on it.
      int status = bind(_server, (struct sockaddr*)_server_addr, sizeof(_server_addr));
      if (status >= 0)
      {
        // Fire up sodium and register our shutdown function.
        if (sodium_init() != -1)
        {
          polymec_atexit(server_stop);

          // Spawn our listener thread.
#if !defined(__STDC_NO_THREADS__) || __STDC_NO_THREADS__
          pthread_attr_init(&_listener_attr);
          pthread_create(&_listener, &_listener_attr, accept_connections, NULL);
#else
          thrd_create(&_listener, accept_connections, NULL);
#endif
        }
        else if (fail)
          polymec_error("polymec_server_start: Couldn't initialize encryption library.");
        else
          log_info("polymec_server_start: Couldn't initialize encryption library. Won't start server.");
      }
      else if (fail)
        polymec_error("polymec_server_start: Couldn't bind port %d for the compute engine server.", port);
      else
        log_info("polymec_server_start: Couldn't bind port %d for the compute engine server. Won't start server.", port);
    }
    else if (fail)
      polymec_error("polymec_server_start: Couldn't open a socket for the compute engine server.");
    else
      log_info("polymec_server_start: Couldn't open a socket for the compute engine server. Won't start server.");
  }
}

