// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_MODEL_H
#define POLYMEC_LUA_MODEL_H

#include "core/lua_types.h"
#include "model/model.h"

// This file contains functions for manipulating model objects in the
// Lua interpreter.

/// \addtogroup model model
///@{

/// \addtogroup lua lua
///@{

/// This function registers the model modules within the interpreter L. It
/// should be called before any of these types are accessed within the
/// interpreter.
int lua_register_model_modules(lua_State* L);

/// Pushes a model m onto L's stack.
void lua_push_model(lua_State* L, model_t* m);

/// Returns true if the item at the given index on L's stack is a
/// model, false if not.
bool lua_is_model(lua_State* L, int index);

/// Returns the model at the given index on L's stack, or NULL if
/// the item there is not a model.
model_t* lua_to_model(lua_State* L, int index);

/// Pushes a probe p onto L's stack.
void lua_push_probe(lua_State* L, probe_t* p);

/// Returns true if the item at the given index on L's stack is a
/// probe, false if not.
bool lua_is_probe(lua_State* L, int index);

/// Returns the probe at the given index on L's stack, or NULL if
/// the item there is not a model probe.
probe_t* lua_to_probe(lua_State* L, int index);

///@}

///@}

#endif

