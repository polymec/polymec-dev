// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_GEOMETRY_H
#define POLYMEC_LUA_GEOMETRY_H

#include "core/lua_types.h"
#include "geometry/coord_mapping.h"
#include "geometry/sd_func.h"
#include "geometry/sdt_func.h"
#include "geometry/polygon.h"
#include "geometry/polyhedron.h"
#include "geometry/point_cloud.h"
#include "geometry/unimesh.h"
#include "geometry/blockmesh.h"
#include "geometry/colmesh.h"
#include "geometry/polymesh.h"
#include "geometry/planar_polymesh.h"

// This file contains functions for manipulating geometric data types in the
// Lua interpreter.

/// \addtogroup geometry geometry
///@{

/// \addtogroup lua lua
///@{

/// This function registers the geometry modules within the interpreter L. It
/// should be called before any of these types are accessed within the
/// interpreter.
int lua_register_geometry_modules(lua_State* L);

/// Pushes a coordinate mapping object X onto L's stack.
void lua_push_coord_mapping(lua_State* L, coord_mapping_t* X);

/// Returns true if the item at the given index on L's stack is a
/// coordinate mapping, false if not.
bool lua_is_coord_mapping(lua_State* L, int index);

/// Returns the coordinate mapping at the given index on L's
/// stack, or NULL if the item there is not a coordinate mapping.
coord_mapping_t* lua_to_coord_mapping(lua_State* L, int index);

/// Pushes a signed distance function f onto L's stack.
void lua_push_sd_func(lua_State* L, sd_func_t* f);

/// Returns true if the item at the given index on L's stack is a
/// signed distance function, false if not.
bool lua_is_sd_func(lua_State* L, int index);

/// Returns the signed distance function at the given index on L's
/// stack, or NULL if the item there is not a signed distance function.
sd_func_t* lua_to_sd_func(lua_State* L, int index);

/// Pushes a temporal signed distance function f onto L's stack.
void lua_push_sdt_func(lua_State* L, sdt_func_t* f);

/// Returns true if the item at the given index on L's stack is a
/// temporal signed distance function, false if not.
bool lua_is_sdt_func(lua_State* L, int index);

/// Returns the temporal signed distance function at the given index
/// on L's stack, or NULL if the item there is not a temporal signed
/// distance function.
sdt_func_t* lua_to_sdt_func(lua_State* L, int index);

/// Pushes a polygon p onto L's stack.
void lua_push_polygon(lua_State* L, polygon_t* p);

/// Returns true if the item at the given index on L's stack is a
/// polygon, false if not.
bool lua_is_polygon(lua_State* L, int index);

/// Returns the polygon at the given index on L's stack, or NULL if the
/// item there is not a polygon.
polygon_t* lua_to_polygon(lua_State* L, int index);

/// Pushes a polyhedron p onto L's stack.
void lua_push_polyhedron(lua_State* L, polyhedron_t* p);

/// Returns true if the item at the given index on L's stack is a
/// polyhedron, false if not.
bool lua_is_polyhedron(lua_State* L, int index);

/// Returns the polyhedron at the given index on L's stack, or NULL if the
/// item there is not a polyhedron.
polyhedron_t* lua_to_polyhedron(lua_State* L, int index);

/// Pushes a tagger t onto L's stack.
void lua_push_tagger(lua_State* L, tagger_t* t);

/// Returns true if the item at the given index on L's stack is a tagger,
/// false if not.
bool lua_is_tagger(lua_State* L, int index);

/// Returns the tagger at the given index on L's stack, or NULL
/// if the item there is not a tagger.
tagger_t* lua_to_tagger(lua_State* L, int index);

/// Pushes a point cloud c onto L's stack.
void lua_push_point_cloud(lua_State* L, point_cloud_t* c);

/// Returns true if the item at the given index on L's stack is a point
/// cloud, false if not.
bool lua_is_point_cloud(lua_State* L, int index);

/// Returns the point cloud at the given index on L's stack, or NULL
/// if the item there is not a point cloud.
point_cloud_t* lua_to_point_cloud(lua_State* L, int index);

/// Pushes a unimesh m onto L's stack.
void lua_push_unimesh(lua_State* L, unimesh_t* m);

/// Returns true if the item at the given index on L's stack is a unimesh,
/// false if not.
bool lua_is_unimesh(lua_State* L, int index);

/// Returns the unimesh at the given index on L's stack, or NULL
/// if the item there is not a unimesh.
unimesh_t* lua_to_unimesh(lua_State* L, int index);

/// Pushes a blockmesh m onto L's stack.
void lua_push_blockmesh(lua_State* L, blockmesh_t* m);

/// Returns true if the item at the given index on L's stack is a blockmesh,
/// false if not.
bool lua_is_blockmesh(lua_State* L, int index);

/// Returns the blockmesh at the given index on L's stack, or NULL
/// if the item there is not a unimesh.
blockmesh_t* lua_to_blockmesh(lua_State* L, int index);

/// Pushes a colmesh m onto L's stack.
void lua_push_colmesh(lua_State* L, colmesh_t* m);

/// Returns true if the item at the given index on L's stack is a colmesh,
/// false if not.
bool lua_is_colmesh(lua_State* L, int index);

/// Returns the colmesh at the given index on L's stack, or NULL
/// if the item there is not a colmesh.
colmesh_t* lua_to_colmesh(lua_State* L, int index);

/// Pushes a polymesh m onto L's stack.
void lua_push_polymesh(lua_State* L, polymesh_t* m);

/// Returns true if the item at the given index on L's stack is a polymesh,
/// false if not.
bool lua_is_polymesh(lua_State* L, int index);

/// Returns the polymesh at the given index on L's stack, or NULL
/// if the item there is not a polymesh.
polymesh_t* lua_to_polymesh(lua_State* L, int index);

/// Pushes a planar_polymesh m onto L's stack.
void lua_push_planar_polymesh(lua_State* L, planar_polymesh_t* m);

/// Returns true if the item at the given index on L's stack is a
/// planar_polymesh, false if not.
bool lua_is_planar_polymesh(lua_State* L, int index);

/// Returns the planar_polymesh at the given index on L's stack, or NULL
/// if the item there is not a planar_polymesh.
planar_polymesh_t* lua_to_planar_polymesh(lua_State* L, int index);

///@}

///@}

#endif

