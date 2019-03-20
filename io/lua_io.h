// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_IO_H
#define POLYMEC_LUA_IO_H

#include "core/lua_types.h"
#include "io/silo_file.h"

// This file contains functions for manipulating I/O types and functions
// in the Lua interpreter. We attempt to expose these data types in a seamless
// fashion using Lua's prototype-based object-oriented formalism.

/// \addtogroup io io
///@{

/// \addtogroup lua lua
///@{

/// This function registers the I/O modules within the interpreter L. It
/// should be called before any of these types are accessed within the
/// interpreter.
int lua_register_io_modules(lua_State* L);

/// Pushes a SILO file s onto L's stack.
void lua_push_silo_file(lua_State* L, silo_file_t* s);

/// Returns true if the item at the given index on L's stack is a SILO
/// file, false if not.
bool lua_is_silo_file(lua_State* L, int index);

/// Returns the SILO file at the given index on L's stack, or NULL
/// if the item there is not a SILO file.
silo_file_t* lua_to_silo_file(lua_State* L, int index);

///@}

///@}

#endif

