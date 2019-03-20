// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_DRIVER_H
#define POLYMEC_LUA_DRIVER_H

#include <stdbool.h>

// This is a forward declaration of the Lua interpreter state.
typedef struct lua_State lua_State;

/// \addtogroup core core
///@{

/// \addtogroup lua lua
///@{

/// This function implements a main function for a Lua driver.
/// It sets up handlers for SIGINT and SIGTERM to perform a clean shut down.
/// \param [in] argc The number of parameters passed as command line arguments.
/// \param [in] argv An array of command line arguments.
/// \param [in] register_types_and_modules A function that performs any work
///             to register types and modules before starting a Lua interpreter.
/// \returns EXIT_SUCCESS if the interpreter runs to completion without incident,
///          EXIT_FAILURE otherwise.
int lua_driver(int argc,
               char** argv,
               int (*register_types_and_functions)(lua_State* L));

/// This returns true if the current process is running inside lua_driver,
/// false if not.
bool lua_driver_running(void);

/// This returns the full path to the input script being executed by lua_driver,
/// or NULL if there is no such input script (or if lua_driver is not running).
const char* lua_driver_script(void);

///@}

///@}

#endif
