// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
/// register_types_and_modules is a function that calls the above functions to 
/// register types and modules before firing up a Lua interpreter to parse 
/// an input file or operate interactively. Its return value is ignored.
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
