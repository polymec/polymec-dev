// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// cylindrical_1d.c - Implementations of interpreter functions for generating
// 1D cylindrical meshes.

#include <strings.h>
#include "core/mesh.h"
#include "core/interpreter.h"
#include "geometry/symmetry_1d.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int cylindrical_1d_uniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 3)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = cylindrical_1d.uniform(r1, r2, N)");
  }

  // Get the arguments.
  real_t r1 = (real_t)lua_tonumber(lua, 1);
  if (r1 < 0.0)
    return luaL_error(lua, "r1 must be non-negative.");
  real_t r2 = (real_t)lua_tonumber(lua, 2);
  if (r1 >= r2)
    return luaL_error(lua, "r1 must be less than r2.");
  int N = (int)lua_tonumber(lua, 3);
  if (N <= 0)
    return luaL_error(lua, "N must be positive.");

  // Create the mesh.
  mesh_t* mesh = create_uniform_cylindrical_1d_mesh(MPI_COMM_WORLD, r1, r2, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

int cylindrical_1d_logarithmic(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = cylindrical_1d.logarithmic(r1, r2, log_factor, N)");
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
  int N = (int)lua_tonumber(lua, 3);
  if (N <= 0)
    return luaL_error(lua, "N must be positive.");

  // Create the mesh.
  mesh_t* mesh = create_logarithmic_cylindrical_1d_mesh(MPI_COMM_WORLD, r1, r2, log_factor, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

int cylindrical_1d_nonuniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = cylindrical_1d.nonuniform({r1, r2, ..., rN})");
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
  mesh_t* mesh = create_nonuniform_cylindrical_1d_mesh(MPI_COMM_WORLD, rs, N);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);

  return 1;
}

