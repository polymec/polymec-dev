// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// mesh_factory.c - Implementations of interpreter functions for generating
// meshes.

#include <strings.h>
#include "core/polymec.h"
#include "core/array.h"
#include "core/mesh.h"
#include "core/partition_mesh.h"
#include "core/kd_tree.h"
#include "core/unordered_map.h"
#include "core/tuple.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/create_rectilinear_mesh.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int mesh_factory_uniform(lua_State* lua);
int mesh_factory_uniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = mesh_factory.uniform(nx, ny, nz, bounds)");
  }

  // Get the arguments.
  int nx = (int)(lua_tonumber(lua, 1));
  int ny = (int)(lua_tonumber(lua, 2));
  int nz = (int)(lua_tonumber(lua, 3));
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    return luaL_error(lua, "nx, ny, and nz must all be positive.");
  if (!lua_isboundingbox(lua, 4))
    return luaL_error(lua, "bounds must be a bounding box.");

  // Bounding box.
  bbox_t* bbox = lua_toboundingbox(lua, 4);

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the mesh.
  // FIXME: For now, we create the mesh on 1 proc and then distribute.
  // FIXME: This will be changed when we have dynamic repartitioning 
  // FIXME: working.
//  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, bbox);
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, bbox);
  migrator_t* m = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  m = NULL;

  // Tag its faces.
  tag_rectilinear_mesh_faces(mesh, "x1", "x2", "y1", "y2", "z1", "z2");

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

docstring_t* mesh_factory_uniform_doc(void);
docstring_t* mesh_factory_uniform_doc()
{
  return docstring_from_string("mesh_factory.uniform(nx, ny, nz, bounds) - returns a uniform lattice of\n"
                               "  nx x ny x nz cells filling the bounding box bounds. The faces on the\n"
                               "  boundary of the mesh are tagged x1, x2, y1, y2, z1, and z2 as appropriate.");
}

int mesh_factory_rectilinear(lua_State* lua);
int mesh_factory_rectilinear(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_issequence(lua, 1) || !lua_issequence(lua, 2) || !lua_issequence(lua, 3))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = mesh_factory.rectilinear(xs, ys, zs)");
  }

  // Get the arguments.
  int nxs, nys, nzs;
  real_t* xs = lua_tosequence(lua, 1, &nxs);
  real_t* ys = lua_tosequence(lua, 2, &nys);
  real_t* zs = lua_tosequence(lua, 3, &nzs);
  if ((nxs < 2) || (nys < 2) || (nzs < 2))
  {
    return luaL_error(lua, "xs, ys, and zs must all contain at least 2 values.");
  }

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the mesh.
  // FIXME: See above: this will be unnecessary for dynamic repartitioning.
//  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, nxs, ys, nys, zs, nzs);
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_SELF, xs, nxs, ys, nys, zs, nzs);
  migrator_t* m = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  m = NULL;

  // Tag its faces.
  tag_rectilinear_mesh_faces(mesh, "x1", "x2", "y1", "y2", "z1", "z2");

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

docstring_t* mesh_factory_rectilinear_doc(void);
docstring_t* mesh_factory_rectilinear_doc()
{
  return docstring_from_string("mesh_factory.rectilinear(xs, ys, zs) - returns a rectilinear lattice of\n"
                               "  cells whose x, y, and z NODAL coordinates are specified by xs, ys, and\n"
                               "  zs, respectively.");
}

#if 0
int mesh_factory_cubic_lattice_periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
    return luaL_error(lua, "Arguments must be 2 boundary mesh (face) tags.");

  for (int i = 1; i <= 2; ++i)
  {
    if (!lua_isstring(lua, i))
      return luaL_error(lua, "Argument %d must be a face tag.", i);
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);

  // Based on the names of the tags, we can decide which periodic BC 
  // mapping function to use.
  periodic_bc_t* bc = NULL;
  if (!strcmp(tag1, "x1"))
  {
    if (!strcmp(tag2, "x2"))
      return luaL_error(lua, "Periodic boundary maps from 'x1' to '%s' (must be 'x2').", tag2);

    bc = cubic_lattice_x_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "x2"))
  {
    if (!strcmp(tag2, "x1"))
      return luaL_error(lua, "Periodic boundary maps from 'x2' to '%s' (must be 'x1').", tag2);

    bc = cubic_lattice_x_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "y1"))
  {
    if (!strcmp(tag2, "y2"))
      return luaL_error(lua, "Periodic boundary maps from 'y1' to '%s' (must be 'y2').", tag2);

    bc = cubic_lattice_y_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "y2"))
  {
    if (!strcmp(tag2, "y1"))
      return luaL_error(lua, "Periodic boundary maps from 'y2' to '%s' (must be 'y1').", tag2);

    bc = cubic_lattice_y_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "z1"))
  {
    if (!strcmp(tag2, "z2"))
      return luaL_error(lua, "Periodic boundary maps from 'z1' to '%s' (must be 'z2').", tag2);

    bc = cubic_lattice_z_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "z2"))
  {
    if (!strcmp(tag2, "z1"))
      return luaL_error(lua, "Periodic boundary maps from 'z2' to '%s' (must be 'z1').", tag2);

    bc = cubic_lattice_z_periodic_bc_new(tag1, tag2);
  }
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}
#endif

