// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "core/interpreter.h"
#include "poisson/poisson_bc.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int dirichlet_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_isscalarfunction(lua, 1) && !lua_isnumber(lua, 1)))
    return luaL_error(lua, "Invalid arguments. Usage:\nbc = dirichlet_bc(F)\nwhere F is a number or function.");

  // Get the argument. 
  st_func_t* F;
  if (lua_isnumber(lua, 1))
  {
    real_t F0 = (real_t)lua_tonumber(lua, 1);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_toscalarfunction(lua, 1);
  }

  // Create a boundary condition object and push it onto the stack.
  poisson_bc_t* bc = poisson_bc_new(1.0, 0.0, F);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

static int neumann_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_isscalarfunction(lua, 1) && !lua_isnumber(lua, 1)))
    return luaL_error(lua, "Invalid arguments. Usage:\nbc = neumann_bc(F)\nwhere F is a number or function.");

  // Get the argument. 
  st_func_t* F;
  if (lua_isnumber(lua, 1))
  {
    real_t F0 = (real_t)lua_tonumber(lua, 1);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_toscalarfunction(lua, 1);
  }

  // Create a boundary condition object and push it onto the stack.
  poisson_bc_t* bc = poisson_bc_new(0.0, 1.0, F);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

static int robin_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || 
      !lua_isnumber(lua, 1) ||
      !lua_isnumber(lua, 2) ||
      (!lua_isscalarfunction(lua, 3) && !lua_isnumber(lua, 3)))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\nbc = robin_bc(alpha, beta, F)\nwhere F is a number or function.");
  }

  // Get the arguments. 
  real_t alpha = (real_t)lua_tonumber(lua, 1);
  real_t beta = (real_t)lua_tonumber(lua, 2);
  st_func_t* F;
  if (lua_isnumber(lua, 3))
  {
    real_t F0 = (real_t)lua_tonumber(lua, 3);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_toscalarfunction(lua, 3);
  }

  // Create a boundary condition object and push it onto the stack.
  poisson_bc_t* bc = poisson_bc_new(alpha, beta, F);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

void interpreter_register_poisson_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "dirichlet_bc", dirichlet_bc);
  interpreter_register_function(interpreter, "neumann_bc", neumann_bc);
  interpreter_register_function(interpreter, "robin_bc", robin_bc);
}

