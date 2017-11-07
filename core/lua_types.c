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

// Places the docstrings table at the top of the stack.
void lua_get_docstrings(lua_State* L);
void lua_get_docstrings(lua_State* L)
{
  lua_getfield(L, LUA_REGISTRYINDEX, "docstrings");
  if (lua_isnil(L, -1))
  {
    lua_pop(L, 1);

    // Create a table with weak keys that associates these keys (objects) with 
    // values (their docstrings).
    lua_newtable(L);
    lua_newtable(L);
    lua_pushliteral(L, "k");
    lua_setfield(L, -2, "__mode");
    lua_setmetatable(L, -2);
    lua_setfield(L, LUA_REGISTRYINDEX, "docstrings");
    lua_getfield(L, LUA_REGISTRYINDEX, "docstrings");
  }

  ASSERT(lua_istable(L, -1));
}

// Documents the object at the given index with the given docstring.
void lua_set_docstring(lua_State* L, int index, const char* docstring);
void lua_set_docstring(lua_State* L, int index, const char* docstring)
{
  ASSERT(!lua_isnil(L, index));

  if (docstring == NULL) 
    return;

  // Convert to an absolute stack index.
  if (index < 0)
    index = lua_gettop(L) + index + 1;

  // Plop the docstrings table on top of the stack.
  lua_get_docstrings(L);

  // Copy the object at the given index to the top of the stack.
  lua_pushnil(L);
  lua_copy(L, index, -1);

  // Put the docstring up top.
  lua_pushstring(L, docstring);

  // Associate the object with its docstring and pop the docstrings table.
  lua_settable(L, -3);
  lua_pop(L, 1);
}

static int lua_open_module(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  ASSERT(lua_islightuserdata(L, -1));
  lua_module_function* funcs = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_doc");
  const char* module_doc = lua_tostring(L, -1);
  lua_pop(L, 2);

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
  f[num_funcs].name = NULL;
  f[num_funcs].func = NULL;

  luaL_checkversion(L);

  // Register the functions with this module.
  lua_createtable(L, 0, 1);
  luaL_setfuncs(L, f, 0);

  // Document this object.
  lua_set_docstring(L, -1, module_doc);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_doc");

  return 1;
}

void lua_register_module(lua_State* L,
                         const char* module_name,
                         const char* module_doc,
                         lua_module_function funcs[]) 
{
  // Load the module by calling lua_open_module with the right globals set.
  lua_pushlightuserdata(L, (void*)funcs);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  lua_pushstring(L, module_doc);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_doc");
  luaL_requiref(L, module_name, lua_open_module, 1);
}

void lua_register_module_function_table(lua_State* L,
                                        const char* module_name,
                                        const char* table_name,
                                        const char* table_doc,
                                        lua_module_function funcs[]) 
{
  lua_getglobal(L, module_name);

  // Create a table containing the functions.
  int num_funcs = 0;
  while (funcs[num_funcs].name != NULL)
    ++num_funcs;
  luaL_Reg f[num_funcs+1];
  for (int i = 0; i < num_funcs; ++i)
  {
    f[i].name = funcs[i].name;
    f[i].func = funcs[i].func;
  }
  f[num_funcs].name = NULL;
  f[num_funcs].func = NULL;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
  luaL_newlib(L, f);
#pragma clang diagnostic pop

  // Register the table in the module.
  lua_setfield(L, -2, table_name);

  // Document the table.
  lua_set_docstring(L, -1, table_doc);

  // Document each of the functions.
  lua_getfield(L, -2, table_name);
  for (int i = 0; i < num_funcs; ++i)
  {
    lua_getfield(L, -1, funcs[i].name);
    lua_set_docstring(L, -1, funcs[i].doc);
    lua_pop(L, 1);
  }
  lua_pop(L, 1);
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
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_doc");
  const char* class_doc = lua_tostring(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_module_function* functions = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  lua_class_method* methods = lua_touserdata(L, -1);
  lua_pop(L, 4);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_doc");
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

    // Document each of the methods.
    for (int i = 0; i < num_methods; ++i)
    {
      lua_getfield(L, -1, methods[i].name);
      lua_set_docstring(L, -1, methods[i].doc);
      lua_pop(L, 1);
    }
  }

  // Create a containing module for this type and register its functions.
  if (functions != NULL)
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

    // Stick the functions into a new class table.
    luaL_checkversion(L);
    lua_createtable(L, 0, 1);
    luaL_setfuncs(L, f, 0);

    // Document each of the functions.
    for (int i = 0; i < num_funcs; ++i)
    {
      lua_getfield(L, -1, functions[i].name);
      lua_set_docstring(L, -1, functions[i].doc);
      lua_pop(L, 1);
    }
  }
  else
    lua_newtable(L);

  // Document this class.
  lua_set_docstring(L, -1, class_doc);

  return 1;
}

// This helper registers the module on top of the stack in the module contained 
// within that module name, if there are dots in the name. If the name has no 
// dots, this function does nothing.
static void register_in_parent_module(lua_State* L, const char* full_module_name)
{
  if (string_contains(full_module_name, "."))
  {
    // Determine the module name.
    char* dot = strstr(full_module_name, ".");
    char* next_dot = strstr(dot+1, ".");
    while (next_dot != NULL)
    {
      dot = next_dot;
      next_dot = strstr(dot, ".");
    }
    size_t len = dot - full_module_name;
    char module[len+1];
    memcpy(module, full_module_name, sizeof(char) * len);
    module[len] = '\0';

    // Fetch the module from the global table.
    lua_getglobal(L, module);
    if (lua_isnil(L, -1))
    {
      luaL_error(L, "Couldn't register %s in module %s.", 
                 &full_module_name[len+1], module);
    }
    else
    {
      lua_pushvalue(L, -2); // Push our module up top.
      lua_setfield(L, -2, &full_module_name[len+1]); // Register in the module.
    }
  }
}

void lua_register_class(lua_State* L,
                        const char* class_name,
                        const char* class_doc,
                        lua_module_function functions[],
                        lua_class_method methods[])
{
  ASSERT(methods != NULL);

  // Load the type into a module by calling lua_open_class.
  // First, though, we need to stash the following into the global registry:
  //   1. class name
  //   2. class doc
  //   3. functions such as constructors, other static functions.
  //   4. methods
  // so that lua_open_class can retrieve them on the far end.
  lua_pushstring(L, class_name);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  lua_pushstring(L, class_doc);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_doc");
  lua_pushlightuserdata(L, (void*)functions);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_pushlightuserdata(L, (void*)methods);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  luaL_requiref(L, class_name, lua_open_class, 1);

  // Modules with dots in the name should be registered in the right 
  // place.
  register_in_parent_module(L, class_name);
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
  lua_class_t* obj = luaL_testudata(L, index, class_name);
  return (obj != NULL) ? obj->context : NULL;
}

void* lua_check_object(lua_State* L,
                       int index,
                       const char* class_name)
{
  void* obj = lua_to_object(L, index, class_name);
  if (obj == NULL)
  {
    luaL_error(L, "Argument must be a %s.", class_name);
    return NULL;
  }
  else 
    return obj;
}

bool lua_is_object(lua_State* L,
                   int index,
                   const char* class_name)
{
  return (luaL_testudata(L, index, class_name) != NULL);
}

//------------------------------------------------------------------------
//                               Records
//------------------------------------------------------------------------

typedef struct 
{
  char* record_type_name;
  lua_record_field* fields;
  lua_record_metamethod* index;

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

  // Do we have a user-defined __index metamethod?
  if (r->index != NULL)
    return r->index->metamethod(L);

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

// Dictionary of user-defined __index metamethods by record name.
static string_ptr_unordered_map_t* lua_record_index_mms = NULL;

// This function is called at exit to destroy the above global objects.
static void destroy_record_index_mms(void)
{
  if (lua_record_index_mms != NULL)
    string_ptr_unordered_map_free(lua_record_index_mms);
}

static int lua_open_record(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_name");
  const char* record_type_name = lua_tostring(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_doc");
  const char* record_type_doc = lua_tostring(L, -1);
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
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_doc");
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

  // Do the same for our user-defined __index metamethods.
  // if we haven't.
  if (lua_record_index_mms == NULL)
  {
    lua_record_index_mms = string_ptr_unordered_map_new();
    polymec_atexit(destroy_record_index_mms);
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
    int num_mm = 0, index_mm = -1, num_index_mm = 0;
    while (metamethods[num_mm].name != NULL)
    {
      // Register any user-defined __index metamethod.
      if (strcmp(metamethods[num_mm].name, "__index") == 0)
      {
        string_ptr_unordered_map_insert_with_k_dtor(lua_record_index_mms, 
                                                    string_dup(record_type_name),
                                                    &metamethods[num_mm],
                                                    string_free);
        index_mm = num_mm;
        num_index_mm = 1;
      }
      ++num_mm;
    }
    num_mm -= num_index_mm;

    // Add the metamethods, skipping over __index if we have it.
    luaL_Reg mm[num_mm+4];
    for (int i = 0; i < num_mm + num_index_mm; ++i)
    {
      int j = ((index_mm >= 0) && (i >= index_mm)) ? i-1 : i;
      mm[j].name = metamethods[i].name;
      mm[j].func = metamethods[i].metamethod;
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
  if (functions != NULL)
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

    // Document each function.
    for (int i = 0; i < num_funcs; ++i)
    {
      lua_getfield(L, -1, functions[i].name);
      lua_set_docstring(L, -1, functions[i].doc);
      lua_pop(L, 1);
    }
  }
  else
    lua_newtable(L);

  // Document the record's module.
  lua_set_docstring(L, -1, record_type_doc);

  return 1;
}

void lua_register_record_type(lua_State* L,
                              const char* record_type_name,
                              const char* record_type_doc,
                              lua_module_function functions[],
                              lua_record_field fields[],
                              lua_record_metamethod metamethods[])
{
  ASSERT(fields != NULL);
  ASSERT(metamethods != NULL);
  // Load the record into a module by calling lua_open_record.
  // First, though, we need to stash the following into the global registry:
  //   1. record type name
  //   2. record type name
  //   3. functions such as constructors, other static functions.
  //   4. fields
  //   5. metamethods
  // so that lua_open_record can retrieve them on the far end.
  lua_pushstring(L, record_type_name);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_name");
  lua_pushstring(L, record_type_doc);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_type_doc");
  lua_pushlightuserdata(L, (void*)functions);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_functions");
  lua_pushlightuserdata(L, (void*)fields);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_fields");
  lua_pushlightuserdata(L, (void*)metamethods);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_record_metamethods");
  luaL_requiref(L, record_type_name, lua_open_record, 1);

  // Modules with dots in the name should be registered in the right 
  // place.
  register_in_parent_module(L, record_type_name);
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
  void** index_p = string_ptr_unordered_map_get(lua_record_index_mms, (char*)record_type_name);
  if (index_p != NULL)
    r->index = (lua_record_metamethod*)(*index_p);
  else
    r->index = NULL;

  // Assign data.
  r->context = context;
  r->dtor = dtor;
}

void* lua_to_record(lua_State* L,
                    int index,
                    const char* record_type_name)
{
  lua_record_t* r = luaL_testudata(L, index, record_type_name);
  return (r != NULL) ? r->context : NULL;
}

void* lua_check_record(lua_State* L,
                       int index,
                       const char* record_type_name)
{
  void* r = lua_to_record(L, index, record_type_name);
  if (r == NULL)
  {
    luaL_error(L, "Argument must be a %s.", record_type_name);
    return NULL;
  }
  else 
    return r;
}

bool lua_is_record(lua_State* L,
                   int index,
                   const char* record_type_name)
{
  return (luaL_testudata(L, index, record_type_name) != NULL);
}
