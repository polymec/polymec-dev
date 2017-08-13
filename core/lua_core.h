// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_CORE_H
#define POLYMEC_LUA_CORE_H

#include "core/polymec.h"
#include "core/lua_types.h"
#include "core/tensor2.h"
#include "core/sp_func.h"
#include "core/st_func.h"
#include "core/mesh.h"
#include "core/point_cloud.h"

// This file contains functions for manipulating basic data types in the 
// Lua interpreter. We attempt to expose these data types in a seamless 
// fashion using Lua's prototype-based object-oriented formalism.

// This function registers the core modules within the interpreter L. It 
// should be called before any of these types are accessed within the 
// interpreter.
int lua_register_core_modules(lua_State* L);

// Pushes a real number z onto L's stack.
void lua_push_real(lua_State* L, real_t x);

// Returns true if the item at the given index on L's stack is a real
// number, false if not.
bool lua_is_real(lua_State* L, int index);

// Returns the real number at the given index on L's stack, or (0, 0)
// if the item is not a real number.
real_t lua_to_real(lua_State* L, int index);

#ifndef __cplusplus
// Pushes a complex number z onto L's stack.
void lua_push_complex(lua_State* L, complex_t z);

// Returns true if the item at the given index on L's stack is a complex
// number, false if not.
bool lua_is_complex(lua_State* L, int index);

// Returns the complex number at the given index on L's stack, or (0, 0)
// if the item is not a complex number.
complex_t lua_to_complex(lua_State* L, int index);
#endif

// Pushes a (3D) point p onto L's stack.
void lua_push_point(lua_State* L, point_t* p);

// Returns true if the item at the given index on L's stack is a point, 
// false if not.
bool lua_is_point(lua_State* L, int index);

// Returns the point at the given index on L's stack, or NULL if the item 
// there is not a point.
point_t* lua_to_point(lua_State* L, int index);

// Pushes a (3D) vector v onto L's stack.
void lua_push_vector(lua_State* L, vector_t* v);

// Returns true if the item at the given index on L's stack is a vector, 
// false if not.
bool lua_is_vector(lua_State* L, int index);

// Returns the vector at the given index on L's stack, or NULL if the item 
// there is not a vector.
vector_t* lua_to_vector(lua_State* L, int index);

// Pushes a (3D) rank-2 tensor t onto L's stack.
void lua_push_tensor2(lua_State* L, tensor2_t* t);

// Returns true if the item at the given index on L's stack is a rank-2 
// tensor, false if not.
bool lua_is_tensor2(lua_State* L, int index);

// Returns the rank-2 tensor at the given index on L's stack, or NULL if 
// the item there is not a rank-2 tensor.
tensor2_t* lua_to_tensor2(lua_State* L, int index);

// Pushes a (3D) symmetric rank-2 tensor t onto L's stack.
void lua_push_sym_tensor2(lua_State* L, sym_tensor2_t* t);

// Returns true if the item at the given index on L's stack is a symmetric 
// rank-2 tensor, false if not.
bool lua_is_sym_tensor2(lua_State* L, int index);

// Returns the symmetric rank-2 tensor at the given index on L's stack, or 
// NULL if the item there is not a symmetric rank-2 tensor.
sym_tensor2_t* lua_to_sym_tensor2(lua_State* L, int index);

// Pushes an MPI communicator onto L's stack.
void lua_push_mpi_comm(lua_State* L, MPI_Comm comm);

// Returns true if the item at the given index on L's stack is an MPI
// communicator, false if not.
bool lua_is_mpi_comm(lua_State* L, int index);

// Returns the MPI communicator at the given index on L's stack, or 
// NULL if the item there is not an MPI communicator.
MPI_Comm lua_to_mpi_comm(lua_State* L, int index);

// This enumerated type describes data stored in an array within a 
// Lua interpreter.
typedef enum
{                        // Array of...
  LUA_ARRAY_BYTE,        // bytes 
  LUA_ARRAY_INT,         // integers 
  LUA_ARRAY_UINT64,      // 64-bit unsigned integers 
  LUA_ARRAY_INT64,       // 64-bit integers 
  LUA_ARRAY_INDEX,       // indices
  LUA_ARRAY_REAL,        // real numbers
  LUA_ARRAY_COMPLEX,     // complex numbers
  LUA_ARRAY_POINT,       // 3D points
  LUA_ARRAY_VECTOR,      // 3D vectors
  LUA_ARRAY_TENSOR2,     // 3D rank-2 tensors
  LUA_ARRAY_SYM_TENSOR2  // 3D symmetric rank-2 tensors
} lua_array_data_t;

// Pushes an array of the given type onto L's stack.
void lua_push_array(lua_State* L, void* array, lua_array_data_t type);

// Returns true if the item at the given index on L's stack is an array 
// of the given type, false if not.
bool lua_is_array(lua_State* L, int index, lua_array_data_t type);

// Returns the array of the given type at the given index on L's stack, or 
// NULL if the item there is not such an array.
void* lua_to_array(lua_State* L, int index, lua_array_data_t type);

// Pushes a multidimensional ndarray of the given rank, shape, and type onto L's stack.
void lua_push_ndarray(lua_State* L, int rank, size_t* shape, 
                      void* array, lua_array_data_t type);

// Returns true if the item at the given index on L's stack is an ndarray 
// of the given type, false if not.
bool lua_is_ndarray(lua_State* L, int index, lua_array_data_t type);

// Returns the ndarray of the given type at the given index on L's stack, or 
// NULL if the item there is not such an ndarray. If this function returns 
// a non-NULL result, rank stores the rank of the array, and shape stores a
// borrowed pointer to the shape array.
void* lua_to_ndarray(lua_State* L, int index, lua_array_data_t type, 
                     int* rank, size_t** shape);

// Pushes a bounding box b onto L's stack.
void lua_push_bbox(lua_State* L, bbox_t* b);

// Returns true if the item at the given index on L's stack is a bounding box,
// false if not.
bool lua_is_bbox(lua_State* L, int index);

// Returns the bounding box at the given index on L's stack, or NULL if the 
// item there is not a bounding box.
bbox_t* lua_to_bbox(lua_State* L, int index);

// Pushes a spatial function f onto L's stack.
void lua_push_sp_func(lua_State* L, sp_func_t* f);

// Returns true if the item at the given index on L's stack is a spatial 
// function, false if not.
bool lua_is_sp_func(lua_State* L, int index);

// Returns the spatial function at the given index on L's stack, or NULL if the 
// item there is not a spatial function.
sp_func_t* lua_to_sp_func(lua_State* L, int index);

// Pushes a space-time function f onto L's stack.
void lua_push_st_func(lua_State* L, st_func_t* f);

// Returns true if the item at the given index on L's stack is a space-time 
// function, false if not.
bool lua_is_st_func(lua_State* L, int index);

// Returns the space-time function at the given index on L's stack, or NULL 
// if the item there is not a space-time function.
st_func_t* lua_to_st_func(lua_State* L, int index);

// Pushes a mesh m onto L's stack.
void lua_push_mesh(lua_State* L, mesh_t* m);

// Returns true if the item at the given index on L's stack is a mesh,
// false if not.
bool lua_is_mesh(lua_State* L, int index);

// Returns the mesh at the given index on L's stack, or NULL 
// if the item there is not a mesh.
mesh_t* lua_to_mesh(lua_State* L, int index);

// Pushes a point cloud c onto L's stack.
void lua_push_point_cloud(lua_State* L, point_cloud_t* c);

// Returns true if the item at the given index on L's stack is a point
// cloud, false if not.
bool lua_is_point_cloud(lua_State* L, int index);

// Returns the point cloud at the given index on L's stack, or NULL 
// if the item there is not a point cloud.
point_cloud_t* lua_to_point_cloud(lua_State* L, int index);

#endif

