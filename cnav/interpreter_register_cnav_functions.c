// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "cnav/interpreter_register_cnav_functions.h"
#include "core/constant_st_func.h"

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
    double F0 = lua_tonumber(lua, 1);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_toscalarfunction(lua, 1);
  }

  // Create a boundary condition object and push it onto the stack.
  cnav_bc_t* bc = cnav_bc_new(1.0, 0.0, F);
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
    double F0 = lua_tonumber(lua, 1);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_toscalarfunction(lua, 1);
  }

  // Create a boundary condition object and push it onto the stack.
  cnav_bc_t* bc = cnav_bc_new(0.0, 1.0, F);
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
  double alpha = lua_tonumber(lua, 1);
  double beta = lua_tonumber(lua, 2);
  st_func_t* F;
  if (lua_isnumber(lua, 3))
  {
    double F0 = lua_tonumber(lua, 3);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_toscalarfunction(lua, 3);
  }

  // Create a boundary condition object and push it onto the stack.
  cnav_bc_t* bc = cnav_bc_new(alpha, beta, F);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

void interpreter_register_cnav_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "dirichlet_bc", dirichlet_bc);
  interpreter_register_function(interpreter, "neumann_bc", neumann_bc);
  interpreter_register_function(interpreter, "robin_bc", robin_bc);
}

