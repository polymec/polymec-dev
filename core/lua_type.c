// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "core/lua_type.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

typedef struct 
{
  // Type and behaviors.
  char* type_name;
  lua_type_attr* attrs;
  lua_type_method* methods;

  // Data.
  void* context;
  void (*dtor)(void*);
} lua_object_t;

// This is the function that gets called when a lua_object is 
// garbage-collected.
static int lua_object_gc(lua_State* L)
{
  lua_object_t* obj = lua_touserdata(L, 1);
  if ((obj->dtor != NULL) && (obj->context != NULL))
    obj->dtor(obj->context);
  return 0;
}

static int lua_object_index(lua_State* L)
{
  // We are given the object and its attribute/method name.
  lua_object_t* obj = lua_touserdata(L, -1);
  const char* attr_or_method = lua_tostring(L, -2);

  // Is it an attribute or a method? Search attributes first.
  int i = 0;
  while ((obj->attrs[i].name != NULL) &&
         (strcmp(obj->attrs[i].name, attr_or_method) != 0)) ++i;
  if (obj->attrs[i].name != NULL) // It's an attribute!
    return obj->attrs[i].getter(L);

  // FIXME: We probably don't need this here.
  // Maybe it's a method.
  i = 0;
  while ((obj->methods[i].name != NULL) &&
         (strcmp(obj->methods[i].name, attr_or_method) != 0)) ++i;
  if (obj->methods[i].name != NULL) // It is!
    return obj->methods[i].method(L);

  // It is neither an attribute nor a method, so there's nothing to do here.
  return luaL_error(L, "%s is neither an attribute nor a method of %s.",
                    attr_or_method, obj->type_name);
}

static int lua_object_newindex(lua_State* L)
{
  // We are given the object and its attribute name.
  lua_object_t* obj = lua_touserdata(L, -1);
  const char* attr = lua_tostring(L, -2);

  // Search the attributes table.
  int i = 0;
  while ((obj->attrs[i].name != NULL) &&
         (strcmp(obj->attrs[i].name, attr) != 0)) ++i;
  if (obj->attrs[i].name != NULL) // Found it!
  {
    if (obj->attrs[i].getter != NULL) // It's writeable!
      return obj->attrs[i].getter(L);
    else
    {
      return luaL_error(L, "%s is a read-only an attribute of %s.",
                        attr, obj->type_name);
    }
  }

  // It's not an attribute after all.
  return luaL_error(L, "%s is not an attribute of %s.",
                    attr, obj->type_name);
}

// Destructor for attribute dictionary entries.
static void destroy_attrs(void* val)
{
  lua_type_attr* attrs = val;
  int i = 0; 
  while (attrs[i].name != NULL)
  {
    string_free((char*)attrs[i].name);
    ++i;
  }
  polymec_free(attrs);
}

// Destructor for method dictionary entries.
static void destroy_methods(void* val)
{
  lua_type_method* methods = val;
  int i = 0; 
  while (methods[i].name != NULL)
  {
    string_free((char*)methods[i].name);
    ++i;
  }
  polymec_free(methods);
}

static int generic_tostring(lua_State* L)
{
  lua_object_t* obj = lua_touserdata(L, -1);
  lua_pushfstring(L, "%s (%p)", obj->type_name, obj->context);
  return 1;
}

// Dictionary of constructors by type.
static string_ptr_unordered_map_t* lua_type_ctors = NULL;

// Dictionary of attributes by type.
static string_ptr_unordered_map_t* lua_type_attrs = NULL;

// Dictionary of methods by type.
static string_ptr_unordered_map_t* lua_type_methods = NULL;

// This function is called at exit to destroy the above global objects.
static void destroy_lua_dicts(void)
{
  if (lua_type_ctors != NULL)
    string_ptr_unordered_map_free(lua_type_ctors);
  if (lua_type_attrs != NULL)
    string_ptr_unordered_map_free(lua_type_attrs);
  if (lua_type_methods != NULL)
    string_ptr_unordered_map_free(lua_type_methods);
}

static int lua_open_type(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_type_type_name");
  const char* type_name = lua_tostring(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_type_ctor");
  int (*ctor)(lua_State* L) = lua_tocfunction(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_type_attr");
  lua_type_attr* attributes = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_type_methods");
  lua_type_method* methods = lua_touserdata(L, -1);
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_type_tostring");
  int (*tostring)(lua_State* L) = lua_tocfunction(L, -1);
  lua_pop(L, 5);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_name");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_ctor");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_attr");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_methods");
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_tostring");

  // First of all, create our global type dictionaries if we haven't 
  // already.
  if (lua_type_ctors == NULL)
  {
    ASSERT(lua_type_attrs == NULL);
    ASSERT(lua_type_methods == NULL);

    lua_type_ctors = string_ptr_unordered_map_new();
    lua_type_attrs = string_ptr_unordered_map_new();
    lua_type_methods = string_ptr_unordered_map_new();
    polymec_atexit(destroy_lua_dicts);
  }

  // Create entries for this type in the global dictionaries.

  // Constructor.
  luaL_Reg* c = polymec_malloc(sizeof(luaL_Reg) * 2);
  c[0].name = string_dup("new");
  c[0].func = ctor;
  c[1].name = NULL;
  c[1].func = NULL;
  string_ptr_unordered_map_insert_with_kv_dtors(lua_type_ctors, 
                                                string_dup(type_name),
                                                c,
                                                string_free,
                                                destroy_methods);

  // Attributes.
  int num_attrs = 0;
  while (attributes[num_attrs].name != NULL)
    ++num_attrs;
  lua_type_attr* a = polymec_malloc(sizeof(lua_type_attr) * (num_attrs+1));
  for (int i = 0; i < num_attrs; ++i)
  {
    a[i].name = string_dup(attributes[i].name);
    a[i].getter = attributes[i].getter;
    a[i].setter = attributes[i].setter;
  }
  a[num_attrs].name = NULL;
  string_ptr_unordered_map_insert_with_kv_dtors(lua_type_attrs, 
                                                string_dup(type_name),
                                                a,
                                                string_free,
                                                destroy_attrs);

  // Methods.
  int num_methods = 0;
  while (methods[num_methods].name != NULL)
    ++num_methods;
  luaL_Reg* m = polymec_malloc(sizeof(luaL_Reg) * (num_methods+4+1));
  for (int i = 0; i < num_methods; ++i)
  {
    m[i].name = string_dup(methods[i].name);
    m[i].func = methods[i].method;
  }
  // Add __index, __newindex, __tostring, __gc methods.
  m[num_methods].name = string_dup("__index");
  m[num_methods].func = lua_object_index;
  m[num_methods+1].name = string_dup("__newindex");
  m[num_methods+1].func = lua_object_newindex;
  m[num_methods+2].name = string_dup("__tostring");
  m[num_methods+2].func = tostring;
  m[num_methods+3].name = string_dup("__gc");
  m[num_methods+3].func = lua_object_gc;
  m[num_methods+4].name = NULL;
  m[num_methods+4].func = NULL;
  string_ptr_unordered_map_insert_with_kv_dtors(lua_type_methods, 
                                                string_dup(type_name),
                                                m,
                                                string_free,
                                                destroy_methods);

  // Create a metatable for this type and populate it.
  luaL_newmetatable(L, type_name);
  luaL_setfuncs(L, m, 0);

  // Create a containing module for this type and register its constructor.
  luaL_checkversion(L);
  lua_createtable(L, 0, 1);
  luaL_setfuncs(L, c, 0);
  return 1;
}

void lua_register_type(lua_State* L,
                       const char* type_name,
                       int (*ctor)(lua_State*),
                       lua_type_attr attributes[],
                       lua_type_method methods[],
                       int (*tostring)(lua_State*))
{
  ASSERT(ctor != NULL);

  // Load the type into a module by calling lua_open_type.
  // First, though, we need to stash the following into the global registry:
  //   1. type_name
  //   2. ctor
  //   3. attributes
  //   4. methods
  //   5. tostring
  // so that lua_open_type can retrieve them on the far end.
  lua_pushstring(L, type_name);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_type_name");
  lua_pushcfunction(L, ctor);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_ctor");
  lua_pushlightuserdata(L, (void*)attributes);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_attr");
  lua_pushlightuserdata(L, (void*)methods);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_methods");
  lua_pushcfunction(L, (tostring != NULL) ? tostring : generic_tostring);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_type_tostring");
  luaL_requiref(L, type_name, lua_open_type, 1);
}

static int lua_open_module(lua_State* L)
{
  // Get stuff out of the registry.
  lua_getfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  ASSERT(lua_islightuserdata(L, -1));
  lua_module_func* funcs = lua_touserdata(L, -1);
  lua_pop(L, 1);

  // Clean up the registry.
  lua_pushnil(L);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");

  // Module functions.
  int num_funcs = 0;
  while (funcs[num_funcs].name != NULL)
    ++num_funcs;
  luaL_Reg* f = polymec_malloc(sizeof(luaL_Reg) * (num_funcs+1));
  for (int i = 0; i < num_funcs; ++i)
  {
    f[i].name = string_dup(funcs[i].name);
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
                         lua_module_func funcs[]) 
{
  // Load the module by calling lua_open_module with the right globals set.
  lua_pushlightuserdata(L, (void*)funcs);
  lua_setfield(L, LUA_REGISTRYINDEX, "lua_open_module_funcs");
  luaL_requiref(L, module_name, lua_open_module, 1);
}

void lua_push_object(lua_State* L,
                     const char* type_name,
                     void* context,
                     void (*dtor)(void* context))
{
  // Is this type registered?
  if (!string_ptr_unordered_map_contains(lua_type_attrs, (char*)type_name))
    luaL_error(L, "%s is not a registered type.", type_name);

  // Allocate the object in the interpreter.
  lua_object_t* obj = lua_newuserdata(L, sizeof(lua_object_t));

  // Fetch the metatable for this type and assign it to obj.
  luaL_getmetatable(L, type_name);
  lua_setmetatable(L, -2);

  // Assign behaviors.
  obj->type_name = string_dup(type_name);
  obj->attrs = *string_ptr_unordered_map_get(lua_type_attrs, (char*)type_name);
  obj->methods = *string_ptr_unordered_map_get(lua_type_methods, (char*)type_name);

  // Assign data.
  obj->context = context;
  obj->dtor = dtor;
}

void* lua_to_object(lua_State* L,
                    int index,
                    const char* type_name)
{
  lua_object_t* obj = luaL_checkudata(L, index, type_name);
  return (obj != NULL) ? obj->context : NULL;
}

bool lua_is_object(lua_State* L,
                   int index,
                   const char* type_name)
{
  return (luaL_checkudata(L, index, type_name) != NULL);
}

