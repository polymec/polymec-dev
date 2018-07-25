// Copyright (c) 2012-2018, Jeffrey N. Johnson
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

// This type represents a field in a Lua module, with a name, 
// a getter, and a setter (if any). A field must have a getter 
// but doesn't need a setter if it is read only.
typedef struct 
{
  const char* name;
  int (*getter)(lua_State* L);
  int (*setter)(lua_State* L);
} lua_module_field;

// This type represents a function for a Lua module.
typedef struct
{
  const char* name;
  int (*func)(lua_State* L);
  const char* doc;
} lua_module_function;

// Registers a set of fields and functions with the interpreter L in the 
// module with the given name.
void lua_register_module(lua_State* L,
                         const char* module_name,
                         const char* module_doc,
                         lua_module_field fields[],
                         lua_module_function functions[]);

// Registers a set of functions in a named table within a module. Useful for 
// factory functions.
void lua_register_module_function_table(lua_State* L,
                                        const char* module_name,
                                        const char* table_name,
                                        const char* table_doc,
                                        lua_module_function funcs[]);

// This type specifies whether a Lua object is owned by the Lua 
// environment or by some C component.
typedef enum
{
  LUA_OWNED_BY_LUA, // owned by Lua environment
  LUA_OWNED_BY_C,   // owned by C environment
} lua_ownership_t;

// This type represents a field in a Lua class, with a name, 
// a getter, and a setter (if any). A field must have a getter 
// but doesn't need a setter if it is read only.
typedef struct 
{
  const char* name;
  int (*getter)(lua_State* L);
  int (*setter)(lua_State* L);
} lua_class_field;

// This type represents a method for a Lua class.
typedef struct
{
  const char* name;
  int (*method)(lua_State* L);
  const char* doc;
} lua_class_method;

// Registers a new Lua class with the interpreter L, giving it a name, 
// a set of static functions and a set of methods. The functions live in 
// a module named after the type. functions may be NULL; methods cannot be.
// The last argument is a C destructor function that frees the data 
// within the context pointer passed to lua_push_object.
void lua_register_class(lua_State* L,
                        const char* class_name,
                        const char* class_doc,
                        lua_module_function functions[],
                        lua_class_field fields[],
                        lua_class_method methods[],
                        void (*c_dtor)(void* context));

// Pushes a new (polymec) Lua object of the given class to the top of the 
// stack in the interpreter L, associating it with a context pointer and a 
// destructor to be called when the object is garbage-collected. 
void lua_push_object(lua_State* L,
                     const char* class_name,
                     void* context);

// Returns true if the object at the given index in the interpreter is of 
// the type identified by the given type name, false if not.
bool lua_is_object(lua_State* L,
                   int index,
                   const char* class_name);

// Returns the context pointer associated with the (polymec) lua object of the 
// given type at the given index in the interpreter L, or NULL if the value at that index is 
// not of that type.
void* lua_to_object(lua_State* L,
                    int index,
                    const char* class_name);

// This is a wrapper around lua_to_object that throws an error if the value at 
// the given index is not an object of the correct class. Analogous to 
// luaL_checkudata.
void* lua_check_object(lua_State* L,
                       int index,
                       const char* class_name);

// This function transfers the ownership of the object with the given type 
// to C or to Lua. By default, Lua owns all objects pushed to the stack.
void lua_transfer_object(lua_State* L, 
                         int index,
                         const char* class_name,
                         lua_ownership_t ownership);

#endif

