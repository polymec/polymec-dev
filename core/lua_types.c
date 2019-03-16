// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

// Container for storing C destructors in Lua.
typedef struct
{
  void (*c_dtor)(void* context);
} lua_c_dtor_t;

// Dictionary of C destructors by object name.
static string_ptr_unordered_map_t* lua_c_dtors = NULL;

static void destroy_c_dtors(void)
{
  if (lua_c_dtors != NULL)
    string_ptr_unordered_map_free(lua_c_dtors);
}

typedef void (*c_dtor_t)(void*);
static c_dtor_t dtor_for_object(const char* type_name)
{
  lua_c_dtor_t* c_dtor = *string_ptr_unordered_map_get(lua_c_dtors, (char*)type_name);
  return c_dtor->c_dtor;
}

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

static int module_index(lua_State* L)
{
  // We are passed the logging table and a field name.
  ASSERT(lua_istable(L, 1));
  const char* field = lua_tostring(L, 2);

  // Fetch the fields from the metatable.
  lua_getmetatable(L, -2);
  lua_getfield(L, -1, "__fields");
  lua_module_field* fields = lua_touserdata(L, -1);
  ASSERT(fields != NULL);

  if (fields != NULL)
  {
    int i = 0;
    while ((fields[i].name != NULL) &&
           (strcmp(fields[i].name, field) != 0)) ++i;
    if (fields[i].name != NULL) // Bingo!
      return fields[i].getter(L);
  }

  // If we didn't find a field, return nil.
  lua_pushnil(L);
  return 1;
}

static int module_newindex(lua_State* L)
{
  // We are passed the logging table, a field name, and a value.
  ASSERT(lua_istable(L, 1));
  const char* field = lua_tostring(L, 2);

  // Fetch the fields from the metatable.
  lua_getmetatable(L, -3);
  lua_getfield(L, -1, "__fields");
  lua_module_field* fields = lua_touserdata(L, -1);
  ASSERT(fields != NULL);

  int i = 0;
  while ((fields[i].name != NULL) &&
         (strcmp(fields[i].name, field) != 0)) ++i;
  if (fields[i].name != NULL) // Found it!
  {
    if (fields[i].setter != NULL) // It's writeable!
    {
      lua_remove(L, 2); // Remove the field name from the stack.
      return fields[i].setter(L);
    }
    else
      return luaL_error(L, "%s is a read-only module field.", field);
  }

  // It's not a field after all.
  return luaL_error(L, "%s is not a module field.", field);
}

static int module_gc(lua_State* L)
{
  // Free the fields.
  lua_getmetatable(L, -1);
  lua_getfield(L, -1, "__fields");
  lua_module_field* fields = lua_touserdata(L, -1);
  ASSERT(fields != NULL);
  int i = 0;
  while (fields[i].name != NULL)
  {
    string_free((char*)fields[i].name);
    ++i;
  }
  return 0;
}

static int lua_open_module(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  ASSERT(lua_islightuserdata(L, -1));
  lua_module_function* funcs = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_fields");
  lua_class_field* fields = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_doc");
  const char* module_doc = lua_tostring(L, -1);
  lua_pop(L, 3);

  // Module functions.
  if (funcs != NULL)
  {
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
  }
  else
    lua_newtable(L); // create an empty module table.

  // If we have fields, create and populate a metatable.
  if (fields != NULL)
  {
    lua_newtable(L);
    lua_pushcfunction(L, module_index);
    lua_setfield(L, -2, "__index");
    lua_pushcfunction(L, module_newindex);
    lua_setfield(L, -2, "__newindex");
    lua_pushcfunction(L, module_gc);
    lua_setfield(L, -2, "__gc");

    // Create a registry of fields witin the metatable.
    int num_fields = 0;
    while (fields[num_fields].name != NULL)
      ++num_fields;
    lua_module_field* f = lua_newuserdata(L, sizeof(lua_module_field) * (num_fields+1));
    for (int i = 0; i < num_fields; ++i)
    {
      f[i].name = string_dup(fields[i].name);
      f[i].getter = fields[i].getter;
      f[i].setter = fields[i].setter;
    }
    f[num_fields].name = NULL;
    lua_setfield(L, -2, "__fields");

    // Add the field names to the metatable so that dir() can find them.
    for (int i = 0; i < num_fields; ++i)
    {
      lua_pushboolean(L, true);
      lua_setfield(L, -2, fields[i].name);
    }

    // Set the metatable for our module table.
    lua_setmetatable(L, -2);
  }

  // Document this module.
  lua_set_docstring(L, -1, module_doc);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_fields");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_doc");

  return 1;
}

void lua_register_module(lua_State* L,
                         const char* module_name,
                         const char* module_doc,
                         lua_module_field fields[],
                         lua_module_function funcs[])
{
  // Load the module by calling lua_open_module with the right globals set.
  lua_pushlightuserdata(L, (void*)funcs);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  lua_pushlightuserdata(L, (void*)fields);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_fields");
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
  ASSERT(funcs != NULL);
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
  char* name;
  lua_class_field* fields;
  lua_class_method* index;
  lua_class_method* newindex;

  void* context;
  void (*dtor)(void*);
} lua_class_t;

// This is the function that gets called when a lua object is
// garbage-collected.
static int lua_class_gc(lua_State* L)
{
  lua_class_t* obj = lua_touserdata(L, 1);
  string_free(obj->name);
  if ((obj->dtor != NULL) && (obj->context != NULL))
    obj->dtor(obj->context);
  return 0;
}

static int lua_class_index(lua_State* L)
{
  // We are given the object and its attribute name.
  lua_class_t* c = lua_touserdata(L, 1);

  // Do we have a user-defined __index metamethod?
  if (c->index != NULL)
    return c->index->method(L);

  // The field must be a name, not a numeric index. Search the fields table.
  const char* field = lua_tostring(L, 2);
  if (c->fields != NULL)
  {
    int i = 0;
    while ((c->fields[i].name != NULL) &&
           (strcmp(c->fields[i].name, field) != 0)) ++i;
    if (c->fields[i].name != NULL) // Bingo!
      return c->fields[i].getter(L);
  }

  // If we didn't find a field, forward the request to our metatable.
  luaL_getmetatable(L, c->name);
  lua_getfield(L, -1, field);
  return 1;
}

static int lua_class_newindex(lua_State* L)
{
  // We are given the object, its field name, and a new field value.
  lua_class_t* c = lua_touserdata(L, 1);

  // Do we have a user-defined __newindex metamethod?
  if (c->newindex != NULL)
    return c->newindex->method(L);

  // The field must be a name, not a numeric index. Search the fields table.
  const char* field = lua_tostring(L, 2);

  if (c->fields != NULL)
  {
    int i = 0;
    while ((c->fields[i].name != NULL) &&
           (strcmp(c->fields[i].name, field) != 0)) ++i;
    if (c->fields[i].name != NULL) // Found it!
    {
      if (c->fields[i].setter != NULL) // It's writeable!
      {
        lua_remove(L, 2); // Remove the field name from the stack.
        return c->fields[i].setter(L);
      }
      else
        return luaL_error(L, "%s is a read-only field of %s.", field, c->name);
    }
  }

  // It's not a field after all.
  return luaL_error(L, "%s is not a field of %s.", field, c->name);
}

// Destructor for class field dictionary entries.
static void destroy_class_fields(void* val)
{
  lua_class_field* fields = val;
  int i = 0;
  while (fields[i].name != NULL)
  {
    string_free((char*)fields[i].name);
    ++i;
  }
  polymec_free(fields);
}

// Dictionary of fields by class name.
static string_ptr_unordered_map_t* lua_class_fields = NULL;

// This function is called at exit to destroy the above global objects.
static void destroy_class_fields_dict(void)
{
  if (lua_class_fields != NULL)
    string_ptr_unordered_map_free(lua_class_fields);
}

// Dictionaries of user-defined __index and __newindex metamethods by class name.
static string_ptr_unordered_map_t* lua_class_index_mms = NULL;
static string_ptr_unordered_map_t* lua_class_newindex_mms = NULL;

// This function is called at exit to destroy the above global objects.
static void destroy_class_index_mms(void)
{
  if (lua_class_index_mms != NULL)
    string_ptr_unordered_map_free(lua_class_index_mms);
  if (lua_class_newindex_mms != NULL)
    string_ptr_unordered_map_free(lua_class_newindex_mms);
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
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_fields");
  lua_class_field* fields = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  lua_class_method* methods = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_class_c_dtor");
  lua_c_dtor_t* c_dtor = lua_touserdata(L, -1);
  lua_pop(L, 6);

  // Create our global dictionary for C destructors if we haven't.
  if (lua_c_dtors == NULL)
  {
    lua_c_dtors = string_ptr_unordered_map_new();
    polymec_atexit(destroy_c_dtors);
  }

  // Register the given C destructor.
  string_ptr_unordered_map_insert_with_kv_dtors(lua_c_dtors,
                                                string_dup(class_name),
                                                c_dtor,
                                                string_free,
                                                polymec_free);

  // Create a global dictionary for class fields if we haven't.
  if (lua_class_fields == NULL)
  {
    lua_class_fields = string_ptr_unordered_map_new();
    polymec_atexit(destroy_class_fields_dict);
  }

  // Do the same for our user-defined __index and __newindex metamethods.
  // if we haven't.
  if (lua_class_index_mms == NULL)
  {
    lua_class_index_mms = string_ptr_unordered_map_new();
    lua_class_newindex_mms = string_ptr_unordered_map_new();
    polymec_atexit(destroy_class_index_mms);
  }

  // Populate the fields entry for this type.
  int num_fields = 0;
  if (fields != NULL)
  {
    while (fields[num_fields].name != NULL)
      ++num_fields;
    lua_class_field* f = polymec_malloc(sizeof(lua_class_field) * (num_fields+1));
    for (int i = 0; i < num_fields; ++i)
    {
      f[i].name = string_dup(fields[i].name);
      f[i].getter = fields[i].getter;
      f[i].setter = fields[i].setter;
    }
    f[num_fields].name = NULL;
    string_ptr_unordered_map_insert_with_kv_dtors(lua_class_fields,
                                                  string_dup(class_name),
                                                  f,
                                                  string_free,
                                                  destroy_class_fields);
  }

  // Create a metatable for this type and populate it with metamethods.
  if (methods != NULL)
  {
    int num_methods = 0, index_mm = -1, newindex_mm = -1, num_index_mm = 0;
    while (methods[num_methods].name != NULL)
    {
      // Register any user-defined __index and __newindex metamethods.
      if (strcmp(methods[num_methods].name, "__index") == 0)
      {
        string_ptr_unordered_map_insert_with_k_dtor(lua_class_index_mms,
                                                    string_dup(class_name),
                                                    &methods[num_methods],
                                                    string_free);
        index_mm = num_methods;
        ++num_index_mm;
      }
      else if (strcmp(methods[num_methods].name, "__newindex") == 0)
      {
        string_ptr_unordered_map_insert_with_k_dtor(lua_class_newindex_mms,
                                                    string_dup(class_name),
                                                    &methods[num_methods],
                                                    string_free);
        newindex_mm = num_methods;
        ++num_index_mm;
      }
      ++num_methods;
    }
    num_methods -= num_index_mm;

    // Add the methods, skipping over __index if we have it.
    luaL_Reg m[num_methods+4];
    for (int i = 0; i < num_methods + num_index_mm; ++i)
    {
      int j = i;
      if ((index_mm >= 0) && (i >= index_mm)) --j;
      if ((newindex_mm >= 0) && (i >= newindex_mm)) --j;
      m[j].name = methods[i].name;
      m[j].func = methods[i].method;
    }

    // Add __index, __newindex, __gc methods.
    m[num_methods].name = "__index";
    m[num_methods].func = lua_class_index;
    m[num_methods+1].name = "__newindex";
    m[num_methods+1].func = lua_class_newindex;
    m[num_methods+2].name = "__gc";
    m[num_methods+2].func = lua_class_gc;
    m[num_methods+3].name = NULL;
    m[num_methods+3].func = NULL;

    // Make the metatable and populate it.
    luaL_newmetatable(L, class_name);
    luaL_setfuncs(L, m, 0);

    // Document each of the methods.
    for (int i = 0; i < num_methods; ++i)
    {
      lua_getfield(L, -1, methods[i].name);
      lua_set_docstring(L, -1, methods[i].doc);
      lua_pop(L, 1);
    }
  }

  // Add the field names to the metatable so that dir() can find them.
  for (int f = 0; f < num_fields; ++f)
  {
    lua_pushboolean(L, true);
    lua_setfield(L, -2, fields[f].name);
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
    f[num_funcs].name = NULL;
    f[num_funcs].func = NULL;

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

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_doc");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_fields");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_c_dtor");

  return 1;
}

// This helper registers the module on top of the stack in the module contained
// within that module name, if there are dots in the name. If the name has no
// dots, this function does nothing.
static void register_in_parent_module(lua_State* L,
                                      const char* full_module_name)
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
      // The parent module sits at index -1, and our module sits at -2.
      lua_pushstring(L, &full_module_name[len+1]);
      lua_pushvalue(L, -3); // Push our module up top.
      lua_rawset(L, -3); // Register in the module directly.
    }
  }
}

void lua_register_class(lua_State* L,
                        const char* class_name,
                        const char* class_doc,
                        lua_module_function functions[],
                        lua_class_field fields[],
                        lua_class_method methods[],
                        void (*c_dtor)(void* context))
{
  // Load the type into a module by calling lua_open_class.
  // First, though, we need to stash the following into the global registry:
  //   1. class name
  //   2. class docstring
  //   3. functions such as constructors, other static functions
  //   4. fields with getters/setters
  //   5. methods
  //   6. a C destructor function
  // so that lua_open_class can retrieve them on the far end.
  lua_pushstring(L, class_name);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_name");
  lua_pushstring(L, class_doc);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_doc");
  lua_pushlightuserdata(L, (void*)functions);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_functions");
  lua_pushlightuserdata(L, (void*)fields);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_fields");
  lua_pushlightuserdata(L, (void*)methods);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_methods");
  lua_c_dtor_t* my_c_dtor = polymec_malloc(sizeof(lua_c_dtor_t));
  my_c_dtor->c_dtor = c_dtor;
  lua_pushlightuserdata(L, (void*)my_c_dtor);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_class_c_dtor");
  luaL_requiref(L, class_name, lua_open_class, 1);

  // Modules with dots in the name should be registered in the right
  // place.
  register_in_parent_module(L, class_name);
}

void lua_push_object(lua_State* L,
                     const char* class_name,
                     void* context)
{
  // Allocate the object in the interpreter.
  lua_class_t* obj = lua_newuserdata(L, sizeof(lua_class_t));

  // Fetch the metatable for this type and assign it to obj.
  luaL_getmetatable(L, class_name);
  lua_setmetatable(L, -2);

  // Assign behaviors.
  obj->name = string_dup(class_name);
  lua_class_field** fields_p = (lua_class_field**)string_ptr_unordered_map_get(lua_class_fields, (char*)class_name);
  if (fields_p != NULL)
    obj->fields = *fields_p;
  else
    obj->fields = NULL;

  void** index_p = string_ptr_unordered_map_get(lua_class_index_mms, (char*)class_name);
  if (index_p != NULL)
    obj->index = (lua_class_method*)(*index_p);
  else
    obj->index = NULL;

  void** newindex_p = string_ptr_unordered_map_get(lua_class_newindex_mms, (char*)class_name);
  if (newindex_p != NULL)
    obj->newindex = (lua_class_method*)(*newindex_p);
  else
    obj->newindex = NULL;

  // Assign data.
  obj->context = context;

  // By default, Lua owns this object, so find the C destructor associated
  // with it.
  obj->dtor = dtor_for_object(class_name);
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

void lua_transfer_object(lua_State* L,
                         int index,
                         const char* class_name,
                         lua_ownership_t ownership)
{
  lua_class_t* obj = luaL_testudata(L, index, class_name);
  ASSERT(obj != NULL);
  if (ownership == LUA_OWNED_BY_C)
    obj->dtor = NULL; // Not our problem!
  else
    obj->dtor = dtor_for_object(class_name);
}

