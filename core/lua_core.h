// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_CORE_H
#define POLYMEC_LUA_CORE_H

#include "core/lua_types.h"
#include "core/point.h"
#include "core/tensor.h"
#include "core/sp_func.h"
#include "core/st_func.h"
#include "core/mesh.h"
#include "core/point_cloud.h"
#include "core/silo_file.h"

// This file contains functions for manipulating basic data types in the 
// Lua interpreter. We attempt to expose these data types in a seamless 
// fashion using Lua's prototype-based object-oriented formalism.

// This function registers the core modules within the interpreter L. It 
// should be called before any of these types are accessed within the 
// interpreter.
int lua_register_core_modules(lua_State* L);

// Pushes a (3D) point p onto L's stack.
void lua_push_point(lua_State* L, point_t* p);

// Returns true if the item at the given index on L's stack is a point, 
// false if not.
bool lua_is_point(lua_State* L, int index);

// Returns the point at the given index on L's stack, or NULL if the item 
// there is not a point.
point_t* lua_to_point(lua_State* L, int index);

// Returns true if the item at the given index on L's stack is or can be 
// interpreted as a list of points. This is true if the item is a table 
// representing a list of either points or 3-tuples.
bool lua_is_point_list(lua_State* L, int index);

// Returns true if and only if the item at the given index on L's stack is 
// a canonical list of points, i.e. a table of points objects.
bool lua_is_canonical_point_list(lua_State* L, int index);

// Converts the item at the given index on L's stack to canonical point list
// form, which is a table representing a list of point objects. If the item 
// at the given index cannot be interpreted as a point list, this function 
// has no effect.
void lua_canonicalize_point_list(lua_State* L, int index);

// Copies data out of a canonical point list at the given index, placing it 
// into the given array of real numbers. If the item is not a canonical point 
// list, this function has no effect. The array must be able to contain 
// a number of reals equal to 3*lua_len(L, index).
void lua_export_point_list(lua_State* L, int index, real_t* array);

// Copies data into a canonical point list at the given index from the given 
// array of real numbers. The number of reals copied out is determined by 
// the number of points in the canonical point list.
void lua_import_point_list(lua_State* L, int index, real_t* array);

// Pushes a (3D) vector v onto L's stack.
void lua_push_vector(lua_State* L, vector_t* v);

// Returns true if the item at the given index on L's stack is a vector, 
// false if not.
bool lua_is_vector(lua_State* L, int index);

// Returns the vector at the given index on L's stack, or NULL if the item 
// there is not a vector.
vector_t* lua_to_vector(lua_State* L, int index);

// Returns true if the item at the given index on L's stack is or can be 
// interpreted as a list of vectors. This is true if the item is a table 
// representing a list of either vectors or 3-tuples.
bool lua_is_vector_list(lua_State* L, int index);

// Returns true if and only if the item at the given index on L's stack is 
// a canonical list of vectors, i.e. a table of vector objects.
bool lua_is_canonical_vector_list(lua_State* L, int index);

// Converts the item at the given index on L's stack to canonical vector list
// form, which is a table representing a list of vector objects. If the item 
// at the given index cannot be interpreted as a vector list, this function 
// has no effect.
void lua_canonicalize_vector_list(lua_State* L, int index);

// Copies data out of a canonical vector list at the given index, placing it 
// into the given array of real numbers. If the item is not a canonical vector 
// list, this function has no effect. The array must be able to contain 
// a number of reals equal to 3*lua_len(L, index).
void lua_export_vector_list(lua_State* L, int index, real_t* array);

// Copies data into a canonical vector list at the given index from the given 
// array of real numbers. The number of reals copied out is determined by 
// the number of vectors in the canonical vector list.
void lua_import_vector_list(lua_State* L, int index, real_t* array);

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

// Pushes a SILO file s onto L's stack.
void lua_push_silo_file(lua_State* L, silo_file_t* s);

// Returns true if the item at the given index on L's stack is a SILO
// file, false if not.
bool lua_is_silo_file(lua_State* L, int index);

// Returns the SILO file at the given index on L's stack, or NULL 
// if the item there is not a SILO file.
silo_file_t* lua_to_silo_file(lua_State* L, int index);

// Pushes a general tensor t onto L's stack.
void lua_push_tensor(lua_State* L, tensor_t* t);

// Returns true if the item at the given index on L's stack can be 
// interpreted as a tensor (i.e., as a table of numbers or tables, or 
// as a number), false if not.
bool lua_is_tensor(lua_State* L, int index);

// Returns the tensor t at the given index on L's stack, or NULL if 
// the item there is not a tensor.
tensor_t* lua_to_tensor(lua_State* L, int index);

#endif

