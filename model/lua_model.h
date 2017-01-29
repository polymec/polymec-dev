// Copyright (c) 2012-2017, Jeffrey N. Johnson
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

// This function registers the model modules within the interpreter L. It 
// should be called before any of these types are accessed within the 
// interpreter.
int lua_register_model_modules(lua_State* L);

// Pushes a model m onto L's stack.
void lua_push_model(lua_State* L, model_t* m);

// Returns true if the item at the given index on L's stack is a 
// model, false if not.
bool lua_is_model(lua_State* L, int index);

// Returns the model at the given index on L's stack, or NULL if 
// the item there is not a model.
model_t* lua_to_model(lua_State* L, int index);

// Pushes a model probe p onto L's stack.
void lua_push_model_probe(lua_State* L, model_probe_t* p);

// Returns true if the item at the given index on L's stack is a 
// model_probe, false if not.
bool lua_is_model_probe(lua_State* L, int index);

// Returns the model_probe at the given index on L's stack, or NULL if 
// the item there is not a model probe.
model_probe_t* lua_to_model_probe(lua_State* L, int index);

// Pushes a model data channel c onto L's stack.
void lua_push_model_data_channel(lua_State* L, model_data_channel_t* c);

// Returns true if the item at the given index on L's stack is a 
// model_data_channel, false if not.
bool lua_is_model_data_channel(lua_State* L, int index);

// Returns the model_data_channel at the given index on L's stack, or 
// NULL if the item there is not a model data channel.
model_data_channel_t* lua_to_model_data_channel(lua_State* L, int index);

// Pushes a model local output o onto L's stack.
void lua_push_model_local_output(lua_State* L, model_local_output_t* o);

// Returns true if the item at the given index on L's stack is a 
// model_local_output, false if not.
bool lua_is_model_local_output(lua_State* L, int index);

// Returns the model_local_output at the given index on L's stack, or 
// NULL if the item there is not a model local output.
model_local_output_t* lua_to_model_local_output(lua_State* L, int index);

#endif

