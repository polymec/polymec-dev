// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// spherical_1d.c - Implementations of interpreter functions for generating
// 1D spherical meshes.

#include <strings.h>
#include "core/mesh.h"
#include "geometry/symmetry.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int spherical_1d_uniform(lua_State* lua);
int spherical_1d_uniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 3)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = spherical_1d.uniform(r1, r2, N)");
  }

  // Get the arguments.
  real_t r1 = (real_t)lua_tonumber(lua, 1);
  if (r1 < 0.0)
    return luaL_error(lua, "r1 must be non-negative.");
  real_t r2 = (real_t)lua_tonumber(lua, 2);
  if (r1 >= r2)
    return luaL_error(lua, "r1 must be less than r2.");
  int N = (int)(lua_tonumber(lua, 3));
  if (N <= 0)
    return luaL_error(lua, "N must be positive.");

  // Create the mesh.
  mesh_t* mesh = create_uniform_spherical_1d_mesh(MPI_COMM_WORLD, r1, r2, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

docstring_t* spherical_1d_uniform_doc(void);
docstring_t* spherical_1d_uniform_doc()
{
  return docstring_from_string("spherical_1d.uniform(r1, r2, N) - returns a uniformly-spaced set of\n"
                               "  N radial cells spanning [r1, r2].");
}

int spherical_1d_logarithmic(lua_State* lua);
int spherical_1d_logarithmic(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = spherical_1d.logarithmic(r1, r2, log_factor, N)");
  }

  // Get the arguments.
  real_t r1 = (real_t)lua_tonumber(lua, 1);
  if (r1 < 0.0)
    return luaL_error(lua, "r1 must be non-negative.");
  real_t r2 = (real_t)lua_tonumber(lua, 2);
  if (r1 >= r2)
    return luaL_error(lua, "r1 must be less than r2.");
  real_t log_factor = (real_t)lua_tonumber(lua, 3);
  if (log_factor <= 0.0)
    return luaL_error(lua, "log factor must be positive.");
  int N = (int)(lua_tonumber(lua, 4));
  if (N <= 0)
    return luaL_error(lua, "N must be positive.");

  // Create the mesh.
  mesh_t* mesh = create_logarithmic_spherical_1d_mesh(MPI_COMM_WORLD, r1, r2, log_factor, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

docstring_t* spherical_1d_logarithmic_doc(void);
docstring_t* spherical_1d_logarithmic_doc()
{
  return docstring_from_string("spherical_1d.logarithmic(r1, r2, log_factor, N) - returns a mesh of\n"
                               "  N radial cells whose spacing is logarithmic, spanning [r1, r2].");
}

int spherical_1d_nonuniform(lua_State* lua);
int spherical_1d_nonuniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = spherical_1d.nonuniform({r1, r2, ..., rN})");
  }

  // Get the arguments.
  int N;
  real_t* rs = lua_tosequence(lua, 1, &N);
  if (rs == NULL)
    return luaL_error(lua, "argument must be a list of coordinates.");
  if (rs[0] < 0.0)
    return luaL_error(lua, "First coordinate must be non-negative.");
  for (int i = 1; i < N; ++i)
  {
    if (rs[i] <= rs[i-1])
      return luaL_error(lua, "coordinates must be monotonically increasing.");
  }

  // Create the mesh.
  mesh_t* mesh = create_nonuniform_spherical_1d_mesh(MPI_COMM_WORLD, rs, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

docstring_t* spherical_1d_nonuniform_doc(void);
docstring_t* spherical_1d_nonuniform_doc()
{
  return docstring_from_string("spherical_1d.nonuniform({r1, r2, ..., rN}) - returns a mesh of\n"
                               "  N cells whose radial coordinates are given, and whose other coordinates\n"
                               "  are zero.");
}

