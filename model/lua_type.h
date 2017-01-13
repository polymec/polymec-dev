// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_TYPE_H
#define POLYMEC_LUA_TYPE_H

#include "core/polymec.h"

// This is a forward declaration of the Lua interpreter state, needed to 
// implement functions in Lua.
typedef struct lua_State lua_State;

// This type represents an attribute in a Lua object, with a name, 
// a getter, and a setter (if any).
typedef struct 
{
  const char* name;
  int (*getter)(lua_State* L);
  int (*setter)(lua_State* L);
} lua_type_attr;

// This type represents a method on a Lua object.
typedef struct
{
  const char* name;
  int (*method)(lua_State* L);
} lua_type_method;

// Registers a new Lua type with the interpreter L, giving it a name, a 
// constructor, and a set of attributes and methods. Optionally, a 
// tostring method can be provided to generate a string for an object of 
// this type.
void lua_register_type(lua_State* L,
                       const char* type_name,
                       int (*ctor)(lua_State*),
                       lua_type_attr attributes[],
                       lua_type_method methods[],
                       int (*tostring)(lua_State*));

// Pushes a new (polymec) Lua object of the given type to the top of the stack 
// in the interpreter L, associating it with a context pointer and a destructor 
// to be called when the object is garbage-collected. 
void lua_push_object(lua_State* L,
                     const char* type_name,
                     void* context,
                     void (*dtor)(void* context));

// Returns the context pointer associated with the (polymec) lua object of the 
// given type at the given index in the interpreter L, or NULL if the value at that index is 
// not of that type.
void* lua_to_object(lua_State* L,
                    int index,
                    const char* type_name);

#endif

