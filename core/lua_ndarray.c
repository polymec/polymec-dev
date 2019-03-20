// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

static int ndarray_get_rank(lua_State* L)
{
  lua_ndarray_t* a = lua_to_object(L, 1, "ndarray");
  lua_pushinteger(L, a->rank);
  return 1;
}

static int ndarray_get_shape(lua_State* L)
{
  lua_ndarray_t* a = lua_to_object(L, 1, "ndarray");
  lua_newtable(L);
  for (int i = 0; i < a->rank; ++i)
  {
    lua_pushinteger(L, (lua_Integer)(a->shape[i]));
    lua_seti(L, -2, (int)i);
  }
  return 1;
}

static lua_class_field ndarray_fields[] = {
  {"rank", ndarray_get_rank, NULL},
  {"shape", ndarray_get_shape, NULL},
  {NULL, NULL, NULL}
};

static void ndarray_dtor(void* context)
{
  lua_ndarray_t* a = context;
  polymec_free(a->shape);
  if (a->owns_data && (a->array != NULL))
    polymec_free(a->array);
  polymec_free(a);
}

static int ndarray_tostring(lua_State* L)
{
  lua_ndarray_t* a = lua_to_object(L, 1, "ndarray");

  char shape_str[24*a->rank];
  sprintf(shape_str, "%d", (int)a->shape[0]);
  for (int i = 1; i < a->rank; ++i)
  {
    char dim_str[24+2];
    sprintf(dim_str, " x %d", (int)a->shape[i]);
    strcat(shape_str, dim_str);
  }

  char type_str[24];
  lua_array_get_type_str(a->type, type_str);

  char str[24 + sizeof(shape_str) + sizeof(type_str)];
  sprintf(str, "ndarray(%s, %s)", shape_str, type_str);

  lua_pushstring(L, str);
  return 1;
}

static lua_class_method ndarray_methods[] = {
  {"__tostring", ndarray_tostring, "String representation metamethod"},
  {NULL, NULL, NULL}
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

static lua_module_function ndarray_funcs[] = {
  {"new", ndarray_new, "ndarray.new(shape, type) -> new n-dimensional array"},
  {NULL, NULL, NULL}
};

void lua_register_ndarray(lua_State* L);
void lua_register_ndarray(lua_State* L)
{
  lua_register_class(L, "ndarray",
                     "An N-dimensional array class.",
                     ndarray_funcs,
                     ndarray_fields,
                     ndarray_methods,
                     ndarray_dtor);
}

void lua_push_ndarray(lua_State* L,
                      int rank,
                      size_t* shape,
                      void* array,
                      lua_array_data_t type)
{
  lua_ndarray_t* a = polymec_malloc(sizeof(lua_ndarray_t));
  a->rank = rank;
  a->shape = polymec_malloc(sizeof(size_t) * rank);
  for (int r = 0; r < rank; ++r)
    a->shape[r] = shape[r];
  a->array = array;
  a->type = type;
  a->owns_data = true;
  lua_push_object(L, "ndarray", a);
}

bool lua_is_ndarray(lua_State* L, int index, lua_array_data_t type)
{
  if (lua_is_object(L, index, "ndarray"))
  {
    lua_ndarray_t* a = lua_to_object(L, index, "ndarray");
    return (a->type == type);
  }
  else
    return false;
}

void* lua_to_ndarray(lua_State* L,
                     int index,
                     lua_array_data_t type,
                     int* rank,
                     size_t** shape)
{
  lua_ndarray_t* a = lua_to_object(L, index, "ndarray");
  if ((a != NULL) && (type == a->type))
  {
    *rank = a->rank;
    *shape = a->shape;
    a->owns_data = false;
    return a->array;
  }
  else
    return NULL;
}

