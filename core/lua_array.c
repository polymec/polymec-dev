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

typedef struct
{
  lua_array_data_t type;
  bool owns_data;
  void* array;
  size_t size;
} lua_array_t;

// This is exposed for use by lua_ndarray.
void lua_array_get_type_str(lua_array_data_t type, char* type_str);
void lua_array_get_type_str(lua_array_data_t type, char* type_str)
{
  if (type == LUA_ARRAY_BYTE)
    strcpy(type_str, "byte");
  else if (type == LUA_ARRAY_INT)
    strcpy(type_str, "int");
  else if (type == LUA_ARRAY_INT64)
    strcpy(type_str, "int64");
  else if (type == LUA_ARRAY_UINT64)
    strcpy(type_str, "uint64");
  else if (type == LUA_ARRAY_INDEX)
    strcpy(type_str, "index");
  else if (type == LUA_ARRAY_REAL)
    strcpy(type_str, "real");
  else if (type == LUA_ARRAY_COMPLEX)
    strcpy(type_str, "complex");
  else if (type == LUA_ARRAY_POINT)
    strcpy(type_str, "point");
  else if (type == LUA_ARRAY_VECTOR)
    strcpy(type_str, "vector");
  else if (type == LUA_ARRAY_TENSOR2)
    strcpy(type_str, "tensor2");
  else // if (type == LUA_ARRAY_SYMTENSOR2)
    strcpy(type_str, "symtensor2");
}

// Also exposed for use by lua_ndarray.
size_t lua_array_elem_size(lua_array_data_t type);
size_t lua_array_elem_size(lua_array_data_t type)
{
  size_t size;
  if (type == LUA_ARRAY_BYTE)
    size = sizeof(uint8_t);
  else if (type == LUA_ARRAY_INT)
    size = sizeof(int);
  else if (type == LUA_ARRAY_INT64)
    size = sizeof(int64_t);
  else if (type == LUA_ARRAY_UINT64)
    size = sizeof(uint64_t);
  else if (type == LUA_ARRAY_INDEX)
    size = sizeof(index_t);
  else if (type == LUA_ARRAY_REAL)
    size = sizeof(real_t);
  else if (type == LUA_ARRAY_COMPLEX)
    size = sizeof(complex_t);
  else if (type == LUA_ARRAY_POINT)
    size = sizeof(point_t);
  else if (type == LUA_ARRAY_VECTOR)
    size = sizeof(vector_t);
  else if (type == LUA_ARRAY_TENSOR2)
    size = sizeof(tensor2_t);
  else // if (type == LUA_ARRAY_SYMTENSOR2)
    size = sizeof(symtensor2_t);
  return size;
}

static int array_index(lua_State* L)
{
  lua_array_t* a = lua_to_object(L, 1, "array");
  if (!lua_isinteger(L, 2) && lua_isstring(L, 2) &&
      (strcmp(lua_tostring(L, 2), "type") == 0))
  {
    char type_str[24];
    lua_array_get_type_str(a->type, type_str);
    lua_pushstring(L, type_str);
    return 1;
  }
  else if (!lua_isinteger(L, 2))
    luaL_error(L, "Array index must be an integer.");
  size_t index = (size_t)lua_tointeger(L, 2);
  if ((index < 1) || (index > a->size))
    luaL_error(L, "Invalid index: must be between 1 and %d.", (int)a->size);
  if (a->type == LUA_ARRAY_BYTE)
  {
    byte_array_t* bytes = a->array;
    lua_pushinteger(L, (lua_Integer)(bytes->data[index-1]));
  }
  else if (a->type == LUA_ARRAY_INT)
  {
    int_array_t* ints = a->array;
    lua_pushinteger(L, (lua_Integer)(ints->data[index-1]));
  }
  else if (a->type == LUA_ARRAY_INT64)
  {
    int64_array_t* ints = a->array;
    lua_pushinteger(L, (lua_Integer)(ints->data[index-1]));
  }
  else if (a->type == LUA_ARRAY_UINT64)
  {
    uint64_array_t* ints = a->array;
    lua_pushinteger(L, (lua_Integer)(ints->data[index-1]));
  }
  else if (a->type == LUA_ARRAY_INDEX)
  {
    index_array_t* indices = a->array;
    lua_pushinteger(L, (lua_Integer)(indices->data[index-1]));
  }
  else if (a->type == LUA_ARRAY_REAL)
  {
    real_array_t* reals = a->array;
    lua_pushnumber(L, (double)(reals->data[index-1]));
  }
  else if (a->type == LUA_ARRAY_COMPLEX)
  {
    complex_array_t* complexes = a->array;
    lua_push_complex(L, complexes->data[index-1]);
  }
  else if (a->type == LUA_ARRAY_POINT)
  {
    point_array_t* points = a->array;
    lua_push_point(L, &points->data[index-1]);
  }
  else if (a->type == LUA_ARRAY_VECTOR)
  {
    vector_array_t* vectors = a->array;
    lua_push_vector(L, &vectors->data[index-1]);
  }
  else if (a->type == LUA_ARRAY_TENSOR2)
  {
    tensor2_array_t* tensors = a->array;
    lua_push_tensor2(L, &tensors->data[index-1]);
  }
  else // if (a->type == LUA_ARRAY_SYMTENSOR2)
  {
    symtensor2_array_t* tensors = a->array;
    lua_push_symtensor2(L, &tensors->data[index-1]);
  }
  return 1;
}

static int array_newindex(lua_State* L)
{
  lua_array_t* a = lua_to_object(L, 1, "array");
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Array index must be an integer.");
  size_t index = (size_t)lua_tointeger(L, 2);
  if ((index < 1) || (index > a->size))
    luaL_error(L, "Invalid index: must be between 1 and %d.", (int)a->size);
  if (a->type == LUA_ARRAY_BYTE)
  {
    if (!lua_isinteger(L, 3))
      luaL_error(L, "Invalid value: must be an integer.");
    byte_array_t* bytes = a->array;
    bytes->data[index-1] = (uint8_t)lua_tointeger(L, 3);
  }
  else if (a->type == LUA_ARRAY_INT)
  {
    if (!lua_isinteger(L, 3))
      luaL_error(L, "Invalid value: must be an integer.");
    int_array_t* ints = a->array;
    ints->data[index-1] = (int)lua_tointeger(L, 3);
  }
  else if (a->type == LUA_ARRAY_INT64)
  {
    if (!lua_isinteger(L, 3))
      luaL_error(L, "Invalid value: must be an integer.");
    int64_array_t* ints = a->array;
    ints->data[index-1] = (int64_t)lua_tointeger(L, 3);
  }
  else if (a->type == LUA_ARRAY_UINT64)
  {
    if (!lua_isinteger(L, 3))
      luaL_error(L, "Invalid value: must be an integer.");
    uint64_array_t* ints = a->array;
    ints->data[index-1] = (uint64_t)lua_tointeger(L, 3);
  }
  else if (a->type == LUA_ARRAY_INDEX)
  {
    if (!lua_isinteger(L, 3))
      luaL_error(L, "Invalid value: must be an integer.");
    index_array_t* indices = a->array;
    indices->data[index-1] = (index_t)lua_tointeger(L, 3);
  }
  else if (a->type == LUA_ARRAY_REAL)
  {
    if (!lua_isnumber(L, 3))
      luaL_error(L, "Invalid value: must be a number.");
    real_array_t* reals = a->array;
    reals->data[index-1] = (real_t)lua_tonumber(L, 3);
  }
  else if (a->type == LUA_ARRAY_COMPLEX)
  {
    if (!lua_is_complex(L, 3))
      luaL_error(L, "Invalid value: must be a complex number.");
    complex_array_t* complexes = a->array;
    complexes->data[index-1] = lua_to_complex(L, 3);
  }
  else if (a->type == LUA_ARRAY_POINT)
  {
    if (!lua_is_point(L, 3))
      luaL_error(L, "Invalid value: must be a point.");
    point_array_t* points = a->array;
    points->data[index-1] = *lua_to_point(L, 3);
  }
  else if (a->type == LUA_ARRAY_VECTOR)
  {
    if (!lua_is_vector(L, 3))
      luaL_error(L, "Invalid value: must be a vector.");
    vector_array_t* vectors = a->array;
    vectors->data[index-1] = *lua_to_vector(L, 3);
  }
  else if (a->type == LUA_ARRAY_TENSOR2)
  {
    if (!lua_is_tensor2(L, 3))
      luaL_error(L, "Invalid value: must be a tensor2.");
    tensor2_array_t* tensors = a->array;
    tensors->data[index-1] = *lua_to_tensor2(L, 3);
  }
  else // if (a->type == LUA_ARRAY_SYMTENSOR2)
  {
    if (!lua_is_symtensor2(L, 3))
      luaL_error(L, "Invalid value: must be a symtensor2.");
    symtensor2_array_t* tensors = a->array;
    tensors->data[index-1] = *lua_to_symtensor2(L, 3);
  }
  return 0;
}

static int array_len(lua_State* L)
{
  lua_array_t* a = lua_to_object(L, 1, "array");
  lua_pushinteger(L, (lua_Integer)a->size);
  return 1;
}

static void* array_from_bytes(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;

  byte_array_t* array = byte_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_isinteger(L, -1))
      byte_array_append(array, (int8_t)(lua_tointeger(L, -1)));
    else
    {
      byte_array_free(array);
      array = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return array;
}

static void* array_from_ints(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;

  int_array_t* array = int_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_isinteger(L, -1))
      int_array_append(array, (int)(lua_tointeger(L, -1)));
    else
    {
      int_array_free(array);
      array = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return array;
}

static void* array_from_int64s(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;

  int64_array_t* array = int64_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_isinteger(L, -1))
      int64_array_append(array, (int64_t)(lua_tointeger(L, -1)));
    else
    {
      int64_array_free(array);
      array = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return array;
}

static void* array_from_uint64s(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;

  uint64_array_t* array = uint64_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_isinteger(L, -1))
      uint64_array_append(array, (uint64_t)(lua_tointeger(L, -1)));
    else
    {
      uint64_array_free(array);
      array = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return array;
}

static void* array_from_indices(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;

  index_array_t* array = index_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_isinteger(L, -1))
      index_array_append(array, (index_t)(lua_tointeger(L, -1)));
    else
    {
      index_array_free(array);
      array = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return array;
}

static void* array_from_reals(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;

  real_array_t* array = real_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_isnumber(L, -1))
      real_array_append(array, (real_t)(lua_tonumber(L, -1)));
    else
    {
      real_array_free(array);
      array = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return array;
}

static void* array_from_complexes(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;
  bool is_complex_list = true;
  complex_array_t* complexes = complex_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    complex_t z;
    bool is_complex = lua_is_complex(L, -1);
    if (is_complex)
      z = lua_to_complex(L, -1);
    else
    {
      bool is_2_tuple = (lua_istable(L, -1) && (lua_rawlen(L, -1) == 2));
      if (!is_2_tuple)
      {
        is_complex_list = false;
        break;
      }
      else
      {
        real_t re, im;
        lua_rawgeti(L, -1, 1);
        if (lua_isnumber(L, -1))
        {
          re = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_complex_list = false;
          break;
        }

        lua_rawgeti(L, -1, 2);
        if (lua_isnumber(L, -1))
        {
          im = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_complex_list = false;
          break;
        }
        z = CMPLX(re, im);
      }
      if (!is_complex_list)
        break;
    }
    complex_array_append(complexes, z);
    lua_pop(L, 1);
  }
  if (!is_complex_list)
  {
    complex_array_free(complexes);
    complexes = NULL;
  }
  return complexes;
}

static void* array_from_points(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;
  bool is_point_list = true;
  point_array_t* points = point_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    point_t p;
    bool is_point = lua_is_point(L, -1);
    if (is_point)
      p = *lua_to_point(L, -1);
    else
    {
      bool is_3_tuple = (lua_istable(L, -1) &&
          (lua_rawlen(L, -1) == 3));
      if (!is_3_tuple)
      {
        is_point_list = false;
        break;
      }
      else
      {

        lua_rawgeti(L, -1, 1);
        if (lua_isnumber(L, -1))
        {
          p.x = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_point_list = false;
          break;
        }

        lua_rawgeti(L, -1, 2);
        if (lua_isnumber(L, -1))
        {
          p.y = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_point_list = false;
          break;
        }

        lua_rawgeti(L, -1, 3);
        if (lua_isnumber(L, -1))
        {
          p.z = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_point_list = false;
          break;
        }
      }
      if (!is_point_list)
        break;
    }
    point_array_append(points, p);
    lua_pop(L, 1);
  }
  if (!is_point_list)
  {
    point_array_free(points);
    points = NULL;
  }
  return points;
}

static void* array_from_vectors(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));
  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;
  bool is_vector_list = true;
  vector_array_t* vecs = vector_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    vector_t v;
    bool is_vec = lua_is_vector(L, -1);
    if (is_vec)
      v = *lua_to_vector(L, -1);
    else
    {
      bool is_3_tuple = (lua_istable(L, -1) &&
          (lua_rawlen(L, -1) == 3));
      if (!is_3_tuple)
      {
        is_vector_list = false;
        break;
      }
      else
      {
        lua_rawgeti(L, -1, 1);
        if (lua_isnumber(L, -1))
        {
          v.x = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_vector_list = false;
          break;
        }

        lua_rawgeti(L, -1, 2);
        if (lua_isnumber(L, -1))
        {
          v.y = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_vector_list = false;
          break;
        }

        lua_rawgeti(L, -1, 3);
        if (lua_isnumber(L, -1))
        {
          v.z = (real_t)lua_tonumber(L, -1);
          lua_pop(L, 1);
        }
        else
        {
          is_vector_list = false;
          break;
        }
      }
      if (!is_vector_list)
        break;
    }
    vector_array_append(vecs, v);
    lua_pop(L, 1);
  }
  if (!is_vector_list)
  {
    vector_array_free(vecs);
    vecs = NULL;
  }
  return vecs;
}

static void* array_from_tensor2s(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;
  tensor2_array_t* tensors = tensor2_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_is_tensor2(L, -1))
      tensor2_array_append(tensors, *lua_to_tensor2(L, -1));
    else
    {
      tensor2_array_free(tensors);
      tensors = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return tensors;
}

static void* array_from_symtensor2s(lua_State* L, int index)
{
  ASSERT(lua_istable(L, index));

  size_t len = lua_rawlen(L, index);
  if (len == 0)
    return NULL;
  symtensor2_array_t* tensors = symtensor2_array_new_with_capacity(len);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, index, (lua_Integer)i);
    if (lua_is_symtensor2(L, -1))
      symtensor2_array_append(tensors, *lua_to_symtensor2(L, -1));
    else
    {
      symtensor2_array_free(tensors);
      tensors = NULL;
      break;
    }
    lua_pop(L, 1);
  }
  return tensors;
}

static void* array_from_type(lua_State* L, int index, lua_array_data_t type)
{
  if (type == LUA_ARRAY_BYTE)
    return array_from_bytes(L, index);
  else if (type == LUA_ARRAY_INT)
    return array_from_ints(L, index);
  else if (type == LUA_ARRAY_INT64)
    return array_from_int64s(L, index);
  else if (type == LUA_ARRAY_UINT64)
    return array_from_uint64s(L, index);
  else if (type == LUA_ARRAY_INDEX)
    return array_from_indices(L, index);
  else if (type == LUA_ARRAY_REAL)
    return array_from_reals(L, index);
  else if (type == LUA_ARRAY_COMPLEX)
    return array_from_complexes(L, index);
  else if (type == LUA_ARRAY_POINT)
    return array_from_points(L, index);
  else if (type == LUA_ARRAY_VECTOR)
    return array_from_vectors(L, index);
  else if (type == LUA_ARRAY_TENSOR2)
    return array_from_tensor2s(L, index);
  else // if (type == LUA_ARRAY_SYMTENSOR2)
    return array_from_symtensor2s(L, index);
}

static int array_concat(lua_State* L)
{
  lua_array_t* a = lua_to_object(L, 1, "array");
  lua_array_data_t type = a->type;
  if (!lua_is_array(L, 2, type) && !lua_istable(L, 2))
    luaL_error(L, "Argument must be an array or a table.");
  lua_array_t* a1 = lua_to_object(L, 2, "array");
  void* ar1 = NULL;
  if (a1 == NULL)
    ar1 = array_from_type(L, 2, type);
  else
    ar1 = a1->array;

  if (ar1 == NULL)
  {
    char type_str[24];
    lua_array_get_type_str(type, type_str);
    luaL_error(L, "Argument could not be converted to an array (%s).", type_str);
  }

  // Create a new array with the contents of ar1 concatenated to a->array.
  if (type == LUA_ARRAY_BYTE)
  {
    byte_array_t* d1 = a->array;
    byte_array_t* d2 = ar1;
    byte_array_t* concat = byte_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      byte_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      byte_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_INT)
  {
    int_array_t* d1 = a->array;
    int_array_t* d2 = ar1;
    int_array_t* concat = int_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      int_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      int_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_INT64)
  {
    int64_array_t* d1 = a->array;
    int64_array_t* d2 = ar1;
    int64_array_t* concat = int64_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      int64_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      int64_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_UINT64)
  {
    uint64_array_t* d1 = a->array;
    uint64_array_t* d2 = ar1;
    uint64_array_t* concat = uint64_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      uint64_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      uint64_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_INDEX)
  {
    index_array_t* d1 = a->array;
    index_array_t* d2 = ar1;
    index_array_t* concat = index_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      index_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      index_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_REAL)
  {
    real_array_t* d1 = a->array;
    real_array_t* d2 = ar1;
    real_array_t* concat = real_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      real_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      real_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_COMPLEX)
  {
    complex_array_t* d1 = a->array;
    complex_array_t* d2 = ar1;
    complex_array_t* concat = complex_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      complex_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      complex_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_POINT)
  {
    point_array_t* d1 = a->array;
    point_array_t* d2 = ar1;
    point_array_t* concat = point_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      point_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      point_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_VECTOR)
  {
    vector_array_t* d1 = a->array;
    vector_array_t* d2 = ar1;
    vector_array_t* concat = vector_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      vector_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      vector_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else if (type == LUA_ARRAY_TENSOR2)
  {
    tensor2_array_t* d1 = a->array;
    tensor2_array_t* d2 = ar1;
    tensor2_array_t* concat = tensor2_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      tensor2_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      tensor2_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  else // if (type == LUA_ARRAY_SYMTENSOR2)
  {
    symtensor2_array_t* d1 = a->array;
    symtensor2_array_t* d2 = ar1;
    symtensor2_array_t* concat = symtensor2_array_new();
    for (size_t i = 0; i < d1->size; ++i)
      symtensor2_array_append(concat, d1->data[i]);
    for (size_t i = 0; i < d2->size; ++i)
      symtensor2_array_append(concat, d2->data[i]);
    lua_push_array(L, concat, type, true);
  }
  return 1;
}

static void array_dtor(void* context)
{
  lua_array_t* a = context;
  if (a->owns_data && (a->array != NULL))
  {
    if (a->type == LUA_ARRAY_BYTE)
      byte_array_free((byte_array_t*)a->array);
    else if (a->type == LUA_ARRAY_INT)
      int_array_free((int_array_t*)a->array);
    else if (a->type == LUA_ARRAY_INT64)
      int64_array_free((int64_array_t*)a->array);
    else if (a->type == LUA_ARRAY_UINT64)
      uint64_array_free((uint64_array_t*)a->array);
    else if (a->type == LUA_ARRAY_INDEX)
      index_array_free((index_array_t*)a->array);
    else if (a->type == LUA_ARRAY_REAL)
      real_array_free((real_array_t*)a->array);
    else if (a->type == LUA_ARRAY_COMPLEX)
      complex_array_free((complex_array_t*)a->array);
    else if (a->type == LUA_ARRAY_POINT)
      point_array_free((point_array_t*)a->array);
    else if (a->type == LUA_ARRAY_VECTOR)
      vector_array_free((vector_array_t*)a->array);
    else if (a->type == LUA_ARRAY_TENSOR2)
      tensor2_array_free((tensor2_array_t*)a->array);
    else // if (a->type == LUA_ARRAY_SYMTENSOR2)
      symtensor2_array_free((symtensor2_array_t*)a->array);
  }
  polymec_free(a);
}

static void write_item(lua_State* L, int index, size_t i, char* str)
{
  lua_array_t* a = lua_to_object(L, index, "array");

  // Get the (i+1)th item in the array.
  lua_pushinteger(L, (lua_Integer)(i+1));
  lua_gettable(L, index);

  // Write out the item.
  if (lua_isinteger(L, -1))
  {
    int val = (int)(lua_tointeger(L, -1));
    if (i < a->size-1)
      snprintf(str, 128, "%d, ", val);
    else
      snprintf(str, 128, "%d", val);
  }
  else if (lua_isnumber(L, -1))
  {
    double val = (double)(lua_tointeger(L, -1));
    if (i < a->size-1)
      snprintf(str, 128, "%g, ", val);
    else
      snprintf(str, 128, "%g", val);
  }
  else // handle general case with __tostring metamethod.
  {
    // Get the __tostring metamethod.
    luaL_getmetafield(L, -1, "__tostring");
    ASSERT(lua_isfunction(L, -1));

    // Push the item itself onto the stack so __tostring calls it.
    lua_pushinteger(L, (lua_Integer)(i+1));
    lua_gettable(L, 1);

    // Call the __tostring metamethod.
    lua_call(L, 1, 1);
    ASSERT(lua_isstring(L, -1));

    // Now write down the thing.
    if (i < a->size-1)
      snprintf(str, 128, "%s, ", lua_tostring(L, -1));
    else
      snprintf(str, 128, "%s", lua_tostring(L, -1));

    lua_pop(L, 1);
  }

  lua_pop(L, 1);
}

static int array_tostring(lua_State* L)
{
  lua_array_t* a = lua_to_object(L, 1, "array");

  char type_str[24];
  lua_array_get_type_str(a->type, type_str);
  char* str = NULL;
  size_t first_leg = a->size, second_leg = 0;
  if (a->size <= 100)
    str = polymec_malloc(sizeof(char) * (32 + a->size * 128));
  else
  {
    first_leg = second_leg = 50;
    str = polymec_malloc(sizeof(char) * (32 + 101 * 128));
  }
  sprintf(str, "array(%s): {", type_str);
  for (size_t i = 0; i < first_leg; ++i)
  {
    char stri[129];
    write_item(L, 1, i, stri);
    strcat(str, stri);
  }
  if (second_leg > 0)
  {
    strcat(str, "..., ");
    for (size_t i = a->size-second_leg; i < a->size; ++i)
    {
      char stri[129];
      write_item(L, 1, i, stri);
      strcat(str, stri);
    }
  }
  strcat(str, "}");

  lua_pushstring(L, str);
  polymec_free(str);
  return 1;
}

static lua_class_method array_methods[] = {
  {"__index", array_index, NULL},
  {"__newindex", array_newindex, NULL},
  {"__len", array_len, NULL},
  {"__concat", array_concat, NULL},
  {"__tostring", array_tostring, NULL},
  {NULL, NULL, NULL}
};

// This is exposed for lua_ndarray, too.
lua_array_data_t lua_array_get_type(lua_State* L, const char* type_str, int index);
lua_array_data_t lua_array_get_type(lua_State* L, const char* type_str, int index)
{
  if (strcmp(type_str, "byte") == 0)
    return LUA_ARRAY_BYTE;
  else if (strcmp(type_str, "int") == 0)
    return LUA_ARRAY_INT;
  else if (strcmp(type_str, "int64") == 0)
    return LUA_ARRAY_INT64;
  else if (strcmp(type_str, "uint64") == 0)
    return LUA_ARRAY_UINT64;
  else if (strcmp(type_str, "index") == 0)
    return LUA_ARRAY_INDEX;
  else if (strcmp(type_str, "real") == 0)
    return LUA_ARRAY_REAL;
  else if (strcmp(type_str, "complex") == 0)
    return LUA_ARRAY_COMPLEX;
  else if (strcmp(type_str, "point") == 0)
    return LUA_ARRAY_POINT;
  else if (strcmp(type_str, "vector") == 0)
    return LUA_ARRAY_VECTOR;
  else if (strcmp(type_str, "tensor2") == 0)
    return LUA_ARRAY_TENSOR2;
  else if (strcmp(type_str, "symtensor2") == 0)
    return LUA_ARRAY_SYMTENSOR2;
  else
  {
    luaL_error(L, "Argument %d must be one of the following:\n"
      "byte, int, int64, uint64, index, real, complex, point, vector, tensor2, symtensor2.", index);
    return LUA_ARRAY_INT;
  }
}

static int array_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    return luaL_error(L, "Arguments must be a table of data and a data type.");
  if (!lua_istable(L, 1))
    return luaL_error(L, "Argument 1 must be a table containing data.");
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument 2 must be a data type.");

  const char* type_str = lua_tostring(L, 2);
  lua_array_data_t type = lua_array_get_type(L, type_str, 2);
  void* array = array_from_type(L, 1, type);
  if (array == NULL)
    return luaL_error(L, "Argument 1 must be a table containing %s data.", type_str);

  lua_push_array(L, array, type, true);
  return 1;
}

static lua_module_function array_funcs[] = {
  {"new", array_new, NULL},
  {NULL, NULL, NULL}
};

void lua_register_array(lua_State* L);
void lua_register_array(lua_State* L)
{
  lua_register_class(L, "array",
                     "A type-savvy array class.",
                     array_funcs,
                     NULL,
                     array_methods,
                     array_dtor);

  // Add data type descriptors to the array table.
  lua_getglobal(L, "array");
  lua_pushstring(L, "byte");
  lua_setfield(L, -2, "byte");
  lua_pushstring(L, "int");
  lua_setfield(L, -2, "int");
  lua_pushstring(L, "int64");
  lua_setfield(L, -2, "int64");
  lua_pushstring(L, "uint64");
  lua_setfield(L, -2, "uint64");
  lua_pushstring(L, "index");
  lua_setfield(L, -2, "index");
  lua_pushstring(L, "real");
  lua_setfield(L, -2, "real");
  lua_pushstring(L, "complex");
  lua_setfield(L, -2, "complex");
  lua_pushstring(L, "point");
  lua_setfield(L, -2, "point");
  lua_pushstring(L, "vector");
  lua_setfield(L, -2, "vector");
  lua_pushstring(L, "tensor2");
  lua_setfield(L, -2, "tensor2");
  lua_pushstring(L, "symtensor2");
  lua_setfield(L, -2, "symtensor2");
}

void lua_push_array(lua_State* L, void* array, lua_array_data_t type, bool free_data)
{
  lua_array_t* a = polymec_malloc(sizeof(lua_array_t));
  if (type == LUA_ARRAY_BYTE)
  {
    byte_array_t* bytes = array;
    a->size = bytes->size;
  }
  else if (type == LUA_ARRAY_INT)
  {
    int_array_t* ints = array;
    a->size = ints->size;
  }
  else if (type == LUA_ARRAY_INT64)
  {
    int64_array_t* ints = array;
    a->size = ints->size;
  }
  else if (type == LUA_ARRAY_UINT64)
  {
    uint64_array_t* ints = array;
    a->size = ints->size;
  }
  else if (type == LUA_ARRAY_INDEX)
  {
    index_array_t* indices = array;
    a->size = indices->size;
  }
  else if (type == LUA_ARRAY_REAL)
  {
    real_array_t* reals = array;
    a->size = reals->size;
  }
  else if (type == LUA_ARRAY_COMPLEX)
  {
    complex_array_t* complexes = array;
    a->size = complexes->size;
  }
  else if (type == LUA_ARRAY_POINT)
  {
    point_array_t* points = array;
    a->size = points->size;
  }
  else if (type == LUA_ARRAY_VECTOR)
  {
    vector_array_t* vectors = array;
    a->size = vectors->size;
  }
  else if (type == LUA_ARRAY_TENSOR2)
  {
    tensor2_array_t* tensors = array;
    a->size = tensors->size;
  }
  else // if (type == LUA_ARRAY_SYMTENSOR2)
  {
    symtensor2_array_t* tensors = array;
    a->size = tensors->size;
  }
  a->array = array;
  a->type = type;
  a->owns_data = free_data;
  lua_push_object(L, "array", a);
}

bool lua_is_array(lua_State* L, int index, lua_array_data_t type)
{
  lua_array_t* a = lua_to_object(L, index, "array");
  return ((a != NULL) && (type == a->type));
}

void* lua_to_array(lua_State* L, int index, lua_array_data_t type)
{
  lua_array_t* a = lua_to_object(L, index, "array");
  if ((a != NULL) && (type == a->type))
    return a->array;
  else
    return NULL;
}

size_t lua_array_size(lua_State* L, int index)
{
  lua_array_t* a = lua_to_object(L, index, "array");
  return a->size;
}

