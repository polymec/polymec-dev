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

// This type represents a function for a Lua module.
typedef struct
{
  const char* name;
  int (*func)(lua_State* L);
} lua_module_function;

// Registers a set of functions with the interpreter L in the module with the 
// given name.
void lua_register_module(lua_State* L,
                         const char* module_name,
                         lua_module_function functions[]);

// This type represents a method for a Lua class.
typedef struct
{
  const char* name;
  int (*method)(lua_State* L);
} lua_class_method;

// Registers a new Lua class with the interpreter L, giving it a name, 
// a set of static functions and a set of methods. The functions live in 
// a module named after the type.
void lua_register_class(lua_State* L,
                        const char* type_name,
                        lua_module_function functions[],
                        lua_class_method methods[]);

// Pushes a new (polymec) Lua object of the given class to the top of the 
// stack in the interpreter L, associating it with a context pointer and a 
// destructor to be called when the object is garbage-collected. 
void lua_push_object(lua_State* L,
                     const char* class_name,
                     void* context,
                     void (*dtor)(void* context));

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

// This type represents a field in a Lua record, with a name, 
// a getter, and a setter (if any). A record must have a getter 
// but may not need a setter if it is read only.
typedef struct 
{
  const char* name;
  int (*getter)(lua_State* L);
  int (*setter)(lua_State* L);
} lua_record_field;

// This type represents a metamethod for a Lua record.
typedef struct
{
  const char* name;
  int (*metamethod)(lua_State* L);
} lua_record_metamethod;

// Registers a new Lua record with the interpreter L, giving it a name,
// a set of static functions and a set of fields. The functions live in 
// a module named after the record. By record, we mean a plain old data type 
// (POD), or a Passive Data Structure (PDS), which has only attributes and 
// no methods or dynamic behavior. Optionally, metamethods may be specified
// to extend the richness of the record in expressions.
void lua_register_record_type(lua_State* L,
                              const char* record_type_name,
                              lua_module_function functions[],
                              lua_record_field fields[],
                              lua_record_metamethod metamethods[]);

// Pushes a new (polymec) Lua record of the given type to the top of the 
// stack in the interpreter L, associating it with a context pointer and a 
// destructor to be called when the object is garbage-collected. 
void lua_push_record(lua_State* L,
                     const char* record_type_name,
                     void* context,
                     void (*dtor)(void* context));

// Returns true if the record at the given index in the interpreter is of 
// the type identified by the given type name, false if not.
bool lua_is_record(lua_State* L,
                   int index,
                   const char* record_type_name);

// Returns the context pointer associated with the (polymec) Lua record of the 
// given type at the given index in the interpreter L, or NULL if the value at that index is 
// not of that type.
void* lua_to_record(lua_State* L,
                    int index,
                    const char* record_type_name);

#endif

