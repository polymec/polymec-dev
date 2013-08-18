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

#include "core/interpreter.h"
#include "core/constant_st_func.h"
#include "core/periodic_bc.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This file contains definitions for functions that are bundled with any 
// interpreter instance.

// Creates a point from a set of coordinates.
static int point(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || 
      !lua_isnumber(lua, 1) || !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3))
  {
    return luaL_error(lua, "Arguments must be x, y, z coordinates.");
  }

  double x = lua_tonumber(lua, 1);
  double y = lua_tonumber(lua, 2);
  double z = lua_tonumber(lua, 3);
  point_t* point = point_new(x, y, z);
  lua_pushpoint(lua, point);
  return 1;
}

// Creates a vector from a set of components.
static int vector(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || 
      !lua_isnumber(lua, 1) || !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3))
  {
    return luaL_error(lua, "Arguments must be x, y, z components.");
  }

  double x = lua_tonumber(lua, 1);
  double y = lua_tonumber(lua, 2);
  double z = lua_tonumber(lua, 3);
  vector_t* vec = vector_new(x, y, z);
  lua_pushvector(lua, vec);
  return 1;
}

// Creates a bounding box from a table.
static int bounding_box(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
  {
    if (!lua_istable(lua, 1))
      return luaL_error(lua, "Argument must be a table containing x1, x2, y1, y2, z1, z2 values.");
  }

  // Look for x1, x2, y1, y2, z1, z2 in the table.
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  const char* entries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  for (int i = 0; i < 6; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1); // Reads name from top, replaces with bounds[name].
    if (!lua_isnumber(lua, -1))
    {
      return luaL_error(lua, "Invalid entry for '%s'.\n"
                        "x1, x2, y1, y2, z1, z2, must all be numbers.", entries[i]);
    }
    switch(i)
    {
      case 0: bbox->x1 = lua_tonumber(lua, -1);
              break;
      case 1: bbox->x2 = lua_tonumber(lua, -1);
              break;
      case 2: bbox->y1 = lua_tonumber(lua, -1);
              break;
      case 3: bbox->y2 = lua_tonumber(lua, -1);
              break;
      case 4: bbox->z1 = lua_tonumber(lua, -1);
              break;
      case 5: bbox->z2 = lua_tonumber(lua, -1);
              break;
      default: break;
    }
    lua_pop(lua, 1); 
  }

  // Push the bounding box onto the stack.
  lua_pushboundingbox(lua, bbox);
  return 1;
}

// Creates a constant function from a number or a 3-tuple.
static int constant_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args == 1) // Scalar-valued constant.
  {
    if (!lua_isnumber(lua, 1))
      return luaL_error(lua, "Argument must be a number.");
  }
  else if (num_args == 3) // Vector-valued constant.
  {
    for (int i = 0; i < 3; ++i)
    {
      if (!lua_isnumber(lua, i+1))
        return luaL_error(lua, "Argument %d must be a number.", i);
    }
  }
  else
    return luaL_error(lua, "Argument must be a 1 or 3 numbers.");

  // Get the arguments.
  double args[3];
  for (int i = 0; i < num_args; ++i)
    args[i] = lua_tonumber(lua, i+1);

  // Push a constant function onto the stack.
  st_func_t* func = constant_st_func_new(num_args, args);
  if (num_args == 1)
    lua_pushscalarfunction(lua, func);
  else
    lua_pushvectorfunction(lua, func);
  return 1;
}

// Creates a vector-valued function from 3 scalar functions.
static int vector_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 3)
    return luaL_error(lua, "Arguments must be 3 scalar functions.");

  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isscalarfunction(lua, i))
      return luaL_error(lua, "Argument %d must be a scalar function.", i);
  }

  st_func_t* functions[3];
  for (int i = 1; i <= 3; ++i)
    functions[i-1] = lua_toscalarfunction(lua, i);
  lua_pushvectorfunction(lua, multicomp_st_func_from_funcs("vector function", functions, 3));
  return 1;
}

// Creates a periodic boundary condition from a pair of tags.
static int periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
    return luaL_error(lua, "Arguments must be 2 boundary mesh (face) tags.");

  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isstring(lua, i))
      return luaL_error(lua, "Argument %d must be a face tag.", i);
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);
  periodic_bc_t* bc = periodic_bc_new(tag1, tag2);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

// Evaluates the gradient of a scalar function.
static int grad(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) && (num_args != 2))
    return luaL_error(lua, "grad: Invalid number of arguments. Must be 2 or 3.");
  
  if ((num_args == 2) && (!lua_isscalarfunction(lua, 1) || 
      (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2))))
  {
    return luaL_error(lua, "grad: Arguments must be (F, x).");
  }
  else if ((num_args == 3) && (!lua_isscalarfunction(lua, 1) || 
      (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2)) ||
      !lua_isnumber(lua, 3)))
  {
    return luaL_error(lua, "grad: Arguments must be (F, x, t).");
  }

  st_func_t* f = lua_toscalarfunction(lua, 1);
  if (!st_func_has_deriv(f, 1))
    return luaL_error(lua, "grad: Scalar function '%s' has no derivative.", st_func_name(f));

  double t = 0.0;
  if (num_args == 3)
    t = lua_tonumber(lua, 3);

  if (lua_ispoint(lua, 2))
  {
    point_t* x = lua_topoint(lua, 2);
    double val[3];
    st_func_eval_deriv(f, 1, x, t, val);
    vector_t* vec = vector_new(val[0], val[1], val[2]);
    lua_pushvector(lua, vec);
  }
  else
  {
    int num_points;
    point_t* x = lua_topointlist(lua, 2, &num_points);
    vector_t* vecs = malloc(sizeof(vector_t) * num_points);
    for (int i = 0; i < num_points; ++i)
    {
      double val[3];
      st_func_eval_deriv(f, 1, &x[i], t, val);
      vecs[i].x = val[0], vecs[i].y = val[1]; vecs[i].z = val[2];
    }
    lua_pushvectorlist(lua, vecs, num_points);
  }

  return 1;
}

void interpreter_register_default_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "point", point);
  interpreter_register_function(interp, "vector", vector);
  interpreter_register_function(interp, "bounding_box", bounding_box);
  interpreter_register_function(interp, "constant_function", constant_function);
  interpreter_register_function(interp, "vector_function", vector_function);
  interpreter_register_function(interp, "periodic_bc", periodic_bc);
  interpreter_register_function(interp, "grad", grad);
}

