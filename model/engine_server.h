// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ENGINE_SERVER_H
#define POLYMEC_ENGINE_SERVER_H

#include "core/polymec.h"

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

// Stops the server after all communications with clients finish.
void engine_server_stop(void);

// Immediately stops the server, interrupting any communications.
void engine_server_halt(void);

#endif

