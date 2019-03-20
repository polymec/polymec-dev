// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_CORE_H
#define POLYMEC_LUA_CORE_H

#include "core/polymec.h"
#include "core/lua_types.h"
#include "core/point2.h"
#include "core/tensor2.h"
#include "core/sp_func.h"
#include "core/st_func.h"
#include "core/rng.h"
#include "core/adj_graph.h"

// This file contains functions for manipulating basic data types in the
// Lua interpreter. We attempt to expose these data types in a seamless
// fashion using Lua's prototype-based object-oriented formalism.

/// \addtogroup core core
///@{

/// \addtogroup lua lua
///@{

/// This function registers the core modules within the interpreter L. It
/// should be called before any of these types are accessed within the
/// interpreter.
int lua_register_core_modules(lua_State* L);

/// Pushes a real number z onto L's stack.
void lua_push_real(lua_State* L, real_t x);

/// Returns true if the item at the given index on L's stack is a real
/// number, false if not.
bool lua_is_real(lua_State* L, int index);

/// Returns the real number at the given index on L's stack, or (0, 0)
/// if the item is not a real number.
real_t lua_to_real(lua_State* L, int index);

#ifndef __cplusplus
/// Pushes a complex number z onto L's stack.
void lua_push_complex(lua_State* L, complex_t z);

/// Returns true if the item at the given index on L's stack is a complex
/// number, false if not.
bool lua_is_complex(lua_State* L, int index);

/// Returns the complex number at the given index on L's stack, or (0, 0)
/// if the item is not a complex number.
complex_t lua_to_complex(lua_State* L, int index);
#endif

/// Pushes a (3D) point p onto L's stack.
void lua_push_point(lua_State* L, point_t* p);

/// Returns true if the item at the given index on L's stack is a (3D) point,
/// false if not.
bool lua_is_point(lua_State* L, int index);

/// Returns the (3D) point at the given index on L's stack, or NULL if the item
/// there is not a (3D) point.
point_t* lua_to_point(lua_State* L, int index);

/// Pushes a 2D point p onto L's stack.
void lua_push_point2(lua_State* L, point2_t* p);

/// Returns true if the item at the given index on L's stack is a 2D point,
/// false if not.
bool lua_is_point2(lua_State* L, int index);

/// Returns the 2D point at the given index on L's stack, or NULL if the item
/// there is not a 2D point.
point2_t* lua_to_point2(lua_State* L, int index);

/// Pushes a (3D) vector v onto L's stack.
void lua_push_vector(lua_State* L, vector_t* v);

/// Returns true if the item at the given index on L's stack is a vector,
/// false if not.
bool lua_is_vector(lua_State* L, int index);

/// Returns the vector at the given index on L's stack, or NULL if the item
/// there is not a vector.
vector_t* lua_to_vector(lua_State* L, int index);

/// Pushes a (3D) rank-2 tensor t onto L's stack.
void lua_push_tensor2(lua_State* L, tensor2_t* t);

/// Returns true if the item at the given index on L's stack is a rank-2
/// tensor, false if not.
bool lua_is_tensor2(lua_State* L, int index);

/// Returns the rank-2 tensor at the given index on L's stack, or NULL if
/// the item there is not a rank-2 tensor.
tensor2_t* lua_to_tensor2(lua_State* L, int index);

/// Pushes a (3D) symmetric rank-2 tensor t onto L's stack.
void lua_push_symtensor2(lua_State* L, symtensor2_t* t);

/// Returns true if the item at the given index on L's stack is a symmetric
/// rank-2 tensor, false if not.
bool lua_is_symtensor2(lua_State* L, int index);

/// Returns the symmetric rank-2 tensor at the given index on L's stack, or
/// NULL if the item there is not a symmetric rank-2 tensor.
symtensor2_t* lua_to_symtensor2(lua_State* L, int index);

/// Pushes an MPI communicator onto L's stack.
void lua_push_mpi_comm(lua_State* L, MPI_Comm comm);

/// Returns true if the item at the given index on L's stack is an MPI
/// communicator, false if not.
bool lua_is_mpi_comm(lua_State* L, int index);

/// Returns the MPI communicator at the given index on L's stack, or
/// NULL if the item there is not an MPI communicator.
MPI_Comm lua_to_mpi_comm(lua_State* L, int index);

/// \enum lua_array_data_t
/// This enumerated type describes data stored in an array within a
/// Lua interpreter.
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
  LUA_ARRAY_SYMTENSOR2   // 3D symmetric rank-2 tensors
} lua_array_data_t;

/// Pushes an array of the given type onto L's stack. If free_data is true,
/// the array's data is destroyed upon collection--otherwise the data is
/// assumed to be owned by the caller or another entity.
void lua_push_array(lua_State* L, void* array, lua_array_data_t type,
                    bool free_data);

/// Returns true if the item at the given index on L's stack is an array
/// of the given type, false if not.
bool lua_is_array(lua_State* L, int index, lua_array_data_t type);

/// Returns the array of the given type at the given index on L's stack, or
/// NULL if the item there is not such an array.
void* lua_to_array(lua_State* L, int index, lua_array_data_t type);

/// Returns the size of the array at the given index, or 0 if the item there
/// is not an array.
size_t lua_array_size(lua_State* L, int index);

/// Pushes a multidimensional ndarray of the given rank, shape, and type onto L's stack.
void lua_push_ndarray(lua_State* L, int rank, size_t* shape,
                      void* array, lua_array_data_t type);

/// Returns true if the item at the given index on L's stack is an ndarray
/// of the given type, false if not.
bool lua_is_ndarray(lua_State* L, int index, lua_array_data_t type);

/// Returns the ndarray of the given type at the given index on L's stack, or
/// NULL if the item there is not such an ndarray. If this function returns
/// a non-NULL result, rank stores the rank of the array, and shape stores a
/// borrowed pointer to the shape array.
void* lua_to_ndarray(lua_State* L, int index, lua_array_data_t type,
                     int* rank, size_t** shape);

/// Pushes a bounding box b onto L's stack.
void lua_push_bbox(lua_State* L, bbox_t* b);

/// Returns true if the item at the given index on L's stack is a bounding box,
/// false if not.
bool lua_is_bbox(lua_State* L, int index);

/// Returns the bounding box at the given index on L's stack, or NULL if the
/// item there is not a bounding box.
bbox_t* lua_to_bbox(lua_State* L, int index);

/// Pushes a spatial function f onto L's stack.
void lua_push_sp_func(lua_State* L, sp_func_t* f);

/// Returns true if the item at the given index on L's stack is a spatial
/// function, false if not.
bool lua_is_sp_func(lua_State* L, int index);

/// Returns the spatial function at the given index on L's stack, or NULL if
/// the item there is not a spatial function.
sp_func_t* lua_to_sp_func(lua_State* L, int index);

/// Constructs a spatial function from the Lua function at the given index on
/// L's stack, or returns NULL if the item there is not a Lua function.
/// \param [in] num_comp The number of components in the function.
sp_func_t* lua_as_sp_func(lua_State* L, int index, int num_comp);

/// Pushes a space-time function f onto L's stack.
void lua_push_st_func(lua_State* L, st_func_t* f);

/// Returns true if the item at the given index on L's stack is a space-time
/// function, false if not.
bool lua_is_st_func(lua_State* L, int index);

/// Returns the space-time function at the given index on L's stack, or NULL
/// if the item there is not a space-time function.
st_func_t* lua_to_st_func(lua_State* L, int index);

/// Constructs a space-time function from the Lua function at the given index
/// on L's stack, or returns NULL if the item there is not a Lua function.
/// \param [in] num_comp The number of components in the function.
st_func_t* lua_as_st_func(lua_State* L, int index, int num_comp);

/// Pushes a random number generator r onto L's stack.
void lua_push_rng(lua_State* L, rng_t* r);

/// Returns true if the item at the given index on L's stack is a random
/// number generator, false if not.
bool lua_is_rng(lua_State* L, int index);

/// Returns the random number generator at the given index on L's stack,
/// or NULL if the item there is not a random number generator.
rng_t* lua_to_rng(lua_State* L, int index);

/// Pushes an adjacency graph g onto L's stack.
void lua_push_adj_graph(lua_State* L, adj_graph_t* g);

/// Returns true if the item at the given index on L's stack is an adjacency
/// graph, false if not.
bool lua_is_adj_graph(lua_State* L, int index);

/// Returns the adjacency graph at the given index on L's stack, or NULL
/// if the item there is not an adjacency graph.
adj_graph_t* lua_to_adj_graph(lua_State* L, int index);

///@}

///@}

#endif

