// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// cartesian_1d.c - Implementations of interpreter functions for generating
// 1D Cartesian meshes.

#include <strings.h>
#include "core/mesh.h"
#include "geometry/symmetry.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int cartesian_1d_uniform(lua_State* lua);
int cartesian_1d_uniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 3)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = cartesian_1d.uniform(x1, x2, N)");
  }

  // Get the arguments.
  real_t x1 = (real_t)lua_tonumber(lua, 1);
  real_t x2 = (real_t)lua_tonumber(lua, 2);
  if (x1 >= x2)
    return luaL_error(lua, "x1 must be less than x2.");
  int N = (int)(lua_tonumber(lua, 3));
  if (N <= 0)
    return luaL_error(lua, "N must be positive.");

  // Create the mesh.
  mesh_t* mesh = create_uniform_cartesian_1d_mesh(MPI_COMM_WORLD, x1, x2, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

docstring_t* cartesian_1d_uniform_doc(void);
docstring_t* cartesian_1d_uniform_doc()
{
  return docstring_from_string("cartesian_1d.uniform(x1, x2, N) - returns a uniformly-spaced set of\n"
                               "  N cells spanning [x1, x2].");
}

int cartesian_1d_logarithmic(lua_State* lua);
int cartesian_1d_logarithmic(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = cartesian_1d.logarithmic(x1, x2, log_factor, N)");
  }

  // Get the arguments.
  real_t x1 = (real_t)lua_tonumber(lua, 1);
  real_t x2 = (real_t)lua_tonumber(lua, 2);
  if (x1 >= x2)
    return luaL_error(lua, "x1 must be less than x2.");
  real_t log_factor = (real_t)lua_tonumber(lua, 3);
  if (log_factor <= 0.0)
    return luaL_error(lua, "log factor must be positive.");
  int N = (int)(lua_tonumber(lua, 4));
  if (N <= 0)
    return luaL_error(lua, "N must be positive.");

  // Create the mesh.
  mesh_t* mesh = create_logarithmic_cartesian_1d_mesh(MPI_COMM_WORLD, x1, x2, log_factor, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

docstring_t* cartesian_1d_logarithmic_doc(void);
docstring_t* cartesian_1d_logarithmic_doc()
{
  return docstring_from_string("cartesian_1d.logarithmic(x1, x2, log_factor, N) - returns a mesh of\n"
                               "  N cells whose spacing is logarithmic, spanning [x1, x2].");
}

int cartesian_1d_nonuniform(lua_State* lua);
int cartesian_1d_nonuniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = cartesian_1d.nonuniform({x1, x2, ..., xN})");
  }

  // Get the arguments.
  int N;
  real_t* xs = lua_tosequence(lua, 1, &N);
  if (xs == NULL)
    return luaL_error(lua, "argument must be a list of coordinates.");

  // Create the mesh.
  mesh_t* mesh = create_nonuniform_cartesian_1d_mesh(MPI_COMM_WORLD, xs, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

docstring_t* cartesian_1d_nonuniform_doc(void);
docstring_t* cartesian_1d_nonuniform_doc()
{
  return docstring_from_string("cartesian_1d.nonuniform({x1, x2, ..., xN}) - returns a mesh of\n"
                               "  N cells whose x coordinates are given, and whose y and z coordinates\n"
                               "  are zero.");
}

