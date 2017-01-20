// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <signal.h>

#include "core/unordered_map.h"
#include "core/lua_types.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

//------------------------------------------------------------------------
//                        Module functions
//------------------------------------------------------------------------

static int lua_open_module(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  ASSERT(lua_islightuserdata(L, -1));
  lua_module_function* funcs = lua_touserdata(L, -1);
  lua_pop(L, 1);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");

  // Module functions.
  int num_funcs = 0;
  while (funcs[num_funcs].name != NULL)
    ++num_funcs;
  luaL_Reg f[num_funcs+1];
  for (int i = 0; i < num_funcs; ++i)
  {
    f[i].name = funcs[i].name;
    f[i].func = funcs[i].func;
  }

  // Register the functions with this module and return it.
  luaL_checkversion(L);
  lua_createtable(L, 0, 1);
  luaL_setfuncs(L, f, 0);
  return 1;
}

void lua_register_module(lua_State* L,
                         const char* module_name,
                         lua_module_function funcs[]) 
{
  // Load the module by calling lua_open_module with the right globals set.
  lua_pushlightuserdata(L, (void*)funcs);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  luaL_requiref(L, module_name, lua_open_module, 1);
}

//------------------------------------------------------------------------
//                               Classes
//------------------------------------------------------------------------

typedef struct 
{
  void* context;
  void (*dtor)(void*);
} lua_class_t;

// This is the function that gets called when a lua object is 
// garbage-collected.
static int lua_class_gc(lua_State* L)
{
  lua_class_t* obj = lua_touserdata(L, 1);
  if ((obj->dtor != NULL) && (obj->context != NULL))
    obj->dtor(obj->context);
  return 0;
}

static int lua_open_class(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  const char* class_name = lua_tostring(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_module_function* functions = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  lua_class_method* methods = lua_touserdata(L, -1);
  lua_pop(L, 3);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");

  // Create a metatable for this type and populate it with methods.
  {
    int num_methods = 0;
    while (methods[num_methods].name != NULL)
      ++num_methods;
    luaL_Reg m[num_methods+2];
    for (int i = 0; i < num_methods; ++i)
    {
      m[i].name = methods[i].name;
      m[i].func = methods[i].method;
    }

    // Add __gc method.
    m[num_methods].name = "__gc";
    m[num_methods].func = lua_class_gc;
    m[num_methods+1].name = NULL;
    m[num_methods+1].func = NULL;

    // Make the metatable.
    luaL_newmetatable(L, class_name);
    luaL_setfuncs(L, m, 0);

    // Set metatable.__index = metatable so that 
    // methods are seamless.
    lua_pushstring(L, "__index");
    lua_pushvalue(L, -2);
    lua_settable(L, -3);
  }

  // Create a containing module for this type and register its functions.
  {
    int num_funcs = 0;
    while (functions[num_funcs].name != NULL)
      ++num_funcs;
    luaL_Reg f[num_funcs + 1];
    for (int i = 0; i < num_funcs; ++i)
    {
      f[i].name = functions[i].name;
      f[i].func = functions[i].func;
    }

    luaL_checkversion(L);
    lua_createtable(L, 0, 1);
    luaL_setfuncs(L, f, 0);
  }

  return 1;
}

void lua_register_class(lua_State* L,
                        const char* class_name,
                        lua_module_function functions[],
                        lua_class_method methods[])
{
  // Load the type into a module by calling lua_open_class.
  // First, though, we need to stash the following into the global registry:
  //   1. class name
  //   2. functions such as constructors, other static functions.
  //   3. methods
  // so that lua_open_class can retrieve them on the far end.
  lua_pushstring(L, class_name);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  lua_pushlightuserdata(L, (void*)functions);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_pushlightuserdata(L, (void*)methods);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  luaL_requiref(L, class_name, lua_open_class, 1);
}

void lua_push_object(lua_State* L,
                     const char* class_name,
                     void* context,
                     void (*dtor)(void* context))
{
  // Allocate the object in the interpreter.
  lua_class_t* obj = lua_newuserdata(L, sizeof(lua_class_t));

  // Fetch the metatable for this type and assign it to obj.
  luaL_getmetatable(L, class_name);
  lua_setmetatable(L, -2);

  // Assign data.
  obj->context = context;
  obj->dtor = dtor;
}

void* lua_to_object(lua_State* L,
                    int index,
                    const char* class_name)
{
  lua_class_t* obj = luaL_checkudata(L, index, class_name);
  return (obj != NULL) ? obj->context : NULL;
}

bool lua_is_object(lua_State* L,
                   int index,
                   const char* class_name)
{
  return (luaL_checkudata(L, index, class_name) != NULL);
}

//------------------------------------------------------------------------
//                               Records
//------------------------------------------------------------------------

typedef struct 
{
  char* record_type_name;
  lua_record_field* fields;

  // Data.
  void* context;
  void (*dtor)(void*);
} lua_record_t;

// This is the function that gets called when a lua_record is 
// garbage-collected.
static int lua_record_gc(lua_State* L)
{
  lua_record_t* r = lua_touserdata(L, 1);
  if ((r->dtor != NULL) && (r->context != NULL))
    r->dtor(r->context);
  return 0;
}

static int lua_record_index(lua_State* L)
{
  // We are given the record and its attribute name.
  lua_record_t* r = lua_touserdata(L, 1);
  const char* field = lua_tostring(L, 2);

  int i = 0;
  while ((r->fields[i].name != NULL) &&
         (strcmp(r->fields[i].name, field) != 0)) ++i;
  if (r->fields[i].name != NULL) // Bingo!
    return r->fields[i].getter(L);

  // Return nil if we didn't find anything.
  lua_pushnil(L);
  return 1;
}

static int lua_record_newindex(lua_State* L)
{
  // We are given the record, its field name, and a new field value.
  lua_record_t* r = lua_touserdata(L, 1);
  const char* field = lua_tostring(L, 2);

  // Search the fields table.
  int i = 0;
  while ((r->fields[i].name != NULL) &&
         (strcmp(r->fields[i].name, field) != 0)) ++i;
  if (r->fields[i].name != NULL) // Found it!
  {
    if (r->fields[i].setter != NULL) // It's writeable!
    {
      lua_remove(L, 2); // Remove the field name from the stack.
      return r->fields[i].setter(L);
    }
    else
    {
      return luaL_error(L, "%s is a read-only field of %s.",
                        field, r->record_type_name);
    }
  }

  // It's not a field after all.
  return luaL_error(L, "%s is not a field of %s.",
                    field, r->record_type_name);
}

// Destructor for field dictionary entries.
static void destroy_fields(void* val)
{
  lua_record_field* fields = val;
  int i = 0; 
  while (fields[i].name != NULL)
  {
    string_free((char*)fields[i].name);
    ++i;
  }
  polymec_free(fields);
}

// Dictionary of fields by record name.
static string_ptr_unordered_map_t* lua_record_fields = NULL;

// This function is called at exit to destroy the above global objects.
static void destroy_record_fields(void)
{
  if (lua_record_fields != NULL)
    string_ptr_unordered_map_free(lua_record_fields);
}

static int lua_open_record(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_name");
  const char* record_type_name = lua_tostring(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_record_functions");
  lua_module_function* functions = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_record_fields");
  lua_record_field* fields = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_record_metamethods");
  lua_record_metamethod* metamethods = lua_touserdata(L, -1);
  lua_pop(L, 4);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_name");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_functions");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_fields");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_metamethods");

  // First of all, create our global dictionary for fields 
  // if we haven't.
  if (lua_record_fields == NULL)
  {
    lua_record_fields = string_ptr_unordered_map_new();
    polymec_atexit(destroy_record_fields);
  }

  // Populate the fields entry for this type.
  int num_fields = 0;
  while (fields[num_fields].name != NULL)
    ++num_fields;
  lua_record_field* f = polymec_malloc(sizeof(lua_record_field) * (num_fields+1));
  for (int i = 0; i < num_fields; ++i)
  {
    f[i].name = string_dup(fields[i].name);
    f[i].getter = fields[i].getter;
    f[i].setter = fields[i].setter;
  }
  f[num_fields].name = NULL;
  string_ptr_unordered_map_insert_with_kv_dtors(lua_record_fields, 
                                                string_dup(record_type_name),
                                                f,
                                                string_free,
                                                destroy_fields);

  // Create a metatable for this type and populate it with metamethods.
  if (metamethods != NULL)
  {
    int num_mm = 0;
    while (metamethods[num_mm].name != NULL)
      ++num_mm;
    luaL_Reg mm[num_mm+4];
    for (int i = 0; i < num_mm; ++i)
    {
      mm[i].name = metamethods[i].name;
      mm[i].func = metamethods[i].metamethod;
    }

    // Add __index, __newindex, __gc methods.
    mm[num_mm].name = "__index";
    mm[num_mm].func = lua_record_index;
    mm[num_mm+1].name = "__newindex";
    mm[num_mm+1].func = lua_record_newindex;
    mm[num_mm+2].name = "__gc";
    mm[num_mm+2].func = lua_record_gc;
    mm[num_mm+3].name = NULL;
    mm[num_mm+3].func = NULL;

    // Make the metatable and populate it.
    luaL_newmetatable(L, record_type_name);
    luaL_setfuncs(L, mm, 0);
  }

  // Create a containing module for this record and register its functions.
  {
    int num_funcs = 0;
    while (functions[num_funcs].name != NULL)
      ++num_funcs;
    luaL_Reg funcs[num_funcs + 1];
    for (int i = 0; i < num_funcs; ++i)
    {
      funcs[i].name = functions[i].name;
      funcs[i].func = functions[i].func;
    }

    luaL_checkversion(L);
    lua_createtable(L, 0, 1);
    luaL_setfuncs(L, funcs, 0);
  }

  return 1;
}

void lua_register_record_type(lua_State* L,
                              const char* record_type_name,
                              lua_module_function functions[],
                              lua_record_field fields[],
                              lua_record_metamethod metamethods[])
{
  // Load the record into a module by calling lua_open_record.
  // First, though, we need to stash the following into the global registry:
  //   1. record type name
  //   2. functions such as constructors, other static functions.
  //   3. fields
  //   4. metamethods
  // so that lua_open_record can retrieve them on the far end.
  lua_pushstring(L, record_type_name);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_name");
  lua_pushlightuserdata(L, (void*)functions);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_functions");
  lua_pushlightuserdata(L, (void*)fields);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_fields");
  lua_pushlightuserdata(L, (void*)metamethods);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_metamethods");
  luaL_requiref(L, record_type_name, lua_open_record, 1);
}

void lua_push_record(lua_State* L,
                     const char* record_type_name,
                     void* context,
                     void (*dtor)(void* context))
{
  // Is this record registered?
  if (!string_ptr_unordered_map_contains(lua_record_fields, (char*)record_type_name))
    luaL_error(L, "%s is not a registered record.", record_type_name);

  // Allocate the record in the interpreter.
  lua_record_t* r = lua_newuserdata(L, sizeof(lua_record_t));

  // Fetch the metatable for this type and assign it to r.
  luaL_getmetatable(L, record_type_name);
  lua_setmetatable(L, -2);

  // Assign behaviors.
  r->record_type_name = string_dup(record_type_name);
  r->fields = *string_ptr_unordered_map_get(lua_record_fields, (char*)record_type_name);

  // Assign data.
  r->context = context;
  r->dtor = dtor;
}

void* lua_to_record(lua_State* L,
                    int index,
                    const char* record_type_name)
{
  lua_record_t* r = luaL_checkudata(L, index, record_type_name);
  return (r != NULL) ? r->context : NULL;
}

bool lua_is_record(lua_State* L,
                   int index,
                   const char* record_type_name)
{
  return (luaL_checkudata(L, index, record_type_name) != NULL);
}
