// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_GEOMETRY_H
#define POLYMEC_LUA_GEOMETRY_H

#include "core/lua_types.h"
#include "geometry/coord_mapping.h"
#include "geometry/sd_func.h"
#include "geometry/sdt_func.h"
#include "geometry/polymesh.h"

// This file contains functions for manipulating geometric data types in the 
// Lua interpreter. 

// This function registers the geometry modules within the interpreter L. It 
// should be called before any of these types are accessed within the 
// interpreter.
int lua_register_geometry_modules(lua_State* L);

// Pushes a coordinate mapping object X onto L's stack.
void lua_push_coord_mapping(lua_State* L, coord_mapping_t* X);

// Returns true if the item at the given index on L's stack is a 
// coordinate mapping, false if not.
bool lua_is_coord_mapping(lua_State* L, int index);

// Returns the coordinate mapping at the given index on L's 
// stack, or NULL if the item there is not a coordinate mapping.
coord_mapping_t* lua_to_coord_mapping(lua_State* L, int index);

// Pushes a signed distance function f onto L's stack.
void lua_push_sd_func(lua_State* L, sd_func_t* f);

// Returns true if the item at the given index on L's stack is a 
// signed distance function, false if not.
bool lua_is_sd_func(lua_State* L, int index);

// Returns the signed distance function at the given index on L's 
// stack, or NULL if the item there is not a signed distance function.
sd_func_t* lua_to_sd_func(lua_State* L, int index);

// Pushes a temporal signed distance function f onto L's stack.
void lua_push_sdt_func(lua_State* L, sdt_func_t* f);

// Returns true if the item at the given index on L's stack is a 
// temporal signed distance function, false if not.
bool lua_is_sdt_func(lua_State* L, int index);

// Returns the temporal signed distance function at the given index 
// on L's stack, or NULL if the item there is not a temporal signed 
// distance function.
sdt_func_t* lua_to_sdt_func(lua_State* L, int index);

// Pushes a polymesh m onto L's stack.
void lua_push_polymesh(lua_State* L, polymesh_t* m);

// Returns true if the item at the given index on L's stack is a polymesh,
// false if not.
bool lua_is_polymesh(lua_State* L, int index);

// Returns the polymesh at the given index on L's stack, or NULL 
// if the item there is not a polymesh.
polymesh_t* lua_to_polymesh(lua_State* L, int index);

#endif

