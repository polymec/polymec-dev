// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ENGINE_SERVER_H
#define POLYMEC_ENGINE_SERVER_H

#include "core/array.h"

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

// Pushes a command to the back of the engine's message queue so that it can be 
// sent to a client. Commands are 128 characters or less. The string is not 
// consumed by the queue. If the engine isn't running, this function has no effect.
void engine_server_push_command(const char* command); 

// Pushes an array of bytes to the back of the engine's message queue so that it 
// can be sent to a client. The data is consumed by the queue. If the engine isn't 
// running, this function has no effect.
void engine_server_push_data(byte_array_t* data);

// Pushes a dynamically-allocated array of real-valued numeric data to the back of 
// the engine's message queue so that it can be sent to a client. The array has a rank 
// (0 for a scalar, 1 or more for a "multi-dimensional" array), a shape (an array with 
// the dimensions of each of the indices), and a pointer to real-valued numbers. The 
// array is consumed by the queue.  If the engine isn't running, this function has no 
// effect.
void engine_server_push_real_array(int rank, size_t* shape, real_t* array);

#ifndef __cplusplus
// Pushes a dynamically-allocated array of complex-valued numeric data to the back 
// of the engine's message queue so that it can be sent to a client. The array has a 
// rank (0 for a scalar, 1 or more for a "multi-dimensional" array), a shape (an array 
// with the dimensions of each of the indices), and a pointer to complex-valued numbers. 
// The array is consumed by the queue.  If the engine isn't running, this function has 
// no effect.
void engine_server_push_complex_array(int rank, size_t* shape, complex_t* array);
#endif

#endif

