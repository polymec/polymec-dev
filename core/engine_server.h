// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ENGINE_SERVER_H
#define POLYMEC_ENGINE_SERVER_H

#include <stdbool.h>

// Polymec allows you to create a server for one or more models to communicate
// model data with a client. This "engine" server allows one connection at a 
// time, and only one server may run for a given simulation.

// Returns true if the engine server is running, false if not.
bool engine_server_is_running(void);

// Starts the server running, using the given directory for managing 
// resources associated with it. If the engine server can't find a key pair
// within this resource directory, it will create a new one.
// Call engine_server_is_running() to make sure this call succeeds. Any 
// problems with starting the server are logged.
void engine_server_start(const char* resource_dir);

// Stops the server after all communications with clients finish. If the 
// engine isn't running, this function has no effect.
void engine_server_stop(void);

// Immediately stops the server, interrupting any communications. If the 
// engine isn't running, this function has no effect.
void engine_server_halt(void);

// Given a plain text password, this returns a newly-allocated string 
// containing a hashed password that can be transmitted across a wire.
char* engine_server_hash_password(const char* plain_text_password);

// Inserts a user into the database of authorized users for this engine, assigning 
// that user a hashed password and a public key. The password is assumed to be
// hashed using the same method as provided by engine_server_hash_password().
// Returns true if the user is successfully inserted, false if not. If the engine 
// isn't running, this function has no effect and returns false.
bool engine_server_insert_user(const char* username, 
                               const char* hashed_password,
                               uint8_t* public_key);

// Deletes a user from the database of authorized users for this engine.
// If the engine isn't running, this function has no effect.
void engine_server_delete_user(const char* username);

// Returns a newly allocated array containing a copy of the engine server's public 
// key. If the engine isn't running, this function returns NULL.
uint8_t* engine_server_public_key(void);

// Given a username and a password hashed using the same method as 
// engine_server_hash_password(), this returns a newly-allocated array containing 
// the user's public key if the user is found and the password is correct, NULL otherwise. 
// If the engine isn't running, this function returns NULL.
uint8_t* engine_server_client_key(const char* username, 
                                  const char* hashed_password);

#endif

