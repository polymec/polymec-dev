// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_TYPE_H
#define POLYMEC_LUA_TYPE_H

// This is a forward declaration of the Lua interpreter state.
typedef struct lua_State lua_State;

// This function implements a main function for a Lua driver. 
// register_types_and_modules is a function that calls the above functions to 
// register types and modules before firing up a Lua interpreter to parse 
// an input file or operate interactively.
int lua_driver(int argc,
               char** argv,
               int (*register_types_and_functions)(lua_State* L));

#endif
