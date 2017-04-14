// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

extern void lua_array_get_type_str(lua_array_data_t type, char* type_str);
extern lua_array_data_t lua_array_get_type(lua_State* L, const char* type_str, int index);
extern size_t lua_array_elem_size(lua_array_data_t type);

typedef struct
{
  int rank;
  size_t* shape;
  void* array;
  lua_array_data_t type;
  bool owns_data;
} lua_ndarray_t;

static int ndarray_gc(lua_State* L)
{
  lua_ndarray_t* a = luaL_checkudata(L, 1, "ndarray");
  polymec_free(a->shape);
  if (a->owns_data && (a->array != NULL))
    polymec_free(a->array);
  return 0;
}

static int ndarray_tostring(lua_State* L)
{
  lua_ndarray_t* a = luaL_checkudata(L, 1, "ndarray");

  char shape_str[24*a->rank];
  sprintf(shape_str, "%d", (int)a->shape[0]);
  for (int i = 1; i < a->rank; ++i)
  {
    char dim_str[24+2];
    sprintf(dim_str, "x %d", (int)a->shape[i]);
    strcat(shape_str, dim_str);
  }

  char type_str[24];
  lua_array_get_type_str(a->type, type_str);

  char str[24 + sizeof(shape_str) + sizeof(type_str)];
  sprintf(str, "ndarray(%s, %s)", shape_str, type_str);

  lua_pushstring(L, str);
  return 1;
}

static luaL_Reg ndarray_mm[] = {
  {"__gc", ndarray_gc},
  {"__tostring", ndarray_tostring},
  {NULL, NULL}
};

static int ndarray_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    return luaL_error(L, "Arguments must be a shape table and a data type.");
  if (!lua_istable(L, 1))
    return luaL_error(L, "Argument 1 must be a shape table.");
  if (!lua_isstring(L, 2)) 
    return luaL_error(L, "Argument 2 must be a data type.");

  // Decipher the shape table.
  int rank = (int)lua_rawlen(L, 1);
  size_t shape[rank];
  for (int i = 1; i <= rank; ++i)
  {
    lua_rawgeti(L, 1, (lua_Integer)i);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Shape table must contain only integers.");
    shape[i-1] = (int)lua_tointeger(L, -1);
  }

  // Decipher the data type.
  const char* type_str = lua_tostring(L, 2);
  lua_array_data_t type = lua_array_get_type(L, type_str, 2);

  // Allocate the array itself, and push it onto the stack.
  size_t size = lua_array_elem_size(type);
  size_t n = 1;
  for (int i = 0; i < rank; ++i)
    n *= shape[i];
  void *array = polymec_malloc(size * n);
  lua_push_ndarray(L, rank, shape, array, type);

  return 1;
}

static luaL_Reg ndarray_funcs[] = {
  {"new", ndarray_new},
  {NULL, NULL}
};

static int lua_open_ndarray(lua_State* L)
{
  // Make the metatable and populate it.
  luaL_newmetatable(L, "ndarray");
  luaL_setfuncs(L, ndarray_mm, 0);

  // Register the containing module for the constructor and 
  // the data type descriptors.
  luaL_checkversion(L);
  lua_createtable(L, 0, 1);
  luaL_setfuncs(L, ndarray_funcs, 0);

  return 1;
}

void lua_register_ndarray(lua_State* L);
void lua_register_ndarray(lua_State* L)
{
  luaL_requiref(L, "ndarray", lua_open_ndarray, 1);
}

void lua_push_ndarray(lua_State* L,
                      int rank,
                      size_t* shape,
                      void* array, 
                      lua_array_data_t type)
{
  lua_ndarray_t* a = lua_newuserdata(L, sizeof(lua_ndarray_t));
  luaL_getmetatable(L, "ndarray");
  lua_setmetatable(L, -2);

  a->rank = rank;
  a->shape = polymec_malloc(sizeof(size_t) * rank);
  for (int r = 0; r < rank; ++r)
    a->shape[r] = shape[r];
  a->array = array;
  a->type = type;
  a->owns_data = false;
}

bool lua_is_ndarray(lua_State* L, int index, lua_array_data_t type)
{
  lua_ndarray_t* a = luaL_testudata(L, index, "ndarray");
  return ((a != NULL) && (type == a->type));
}

void* lua_to_ndarray(lua_State* L, int index, lua_array_data_t type, int* rank, size_t** shape)
{
  lua_ndarray_t* a = luaL_testudata(L, index, "ndarray");
  if ((a != NULL) && (type == a->type))
  {
    *rank = a->rank;
    *shape = a->shape;
    return a;
  }
  else
    return NULL;
}

