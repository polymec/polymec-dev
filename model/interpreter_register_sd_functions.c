// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/sphere_sp_func.h"
#include "geometry/rect_prism_sp_func.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static docstring_t* signed_distance_functions_doc()
{
  return docstring_from_string("signed_distance_functions - A family of signed distance functions that can\n"
                               "  be used to construct geometries based on zero level sets.");
}

static int sphere(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispoint(lua, 1) || !lua_isnumber(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "F = sphere(x, R)");
  }

  // Get the arguments.
  point_t* x = lua_topoint(lua, 1);
  real_t r = (real_t)lua_tonumber(lua, 2);
  if (r <= 0.0)
    return luaL_error(lua, "Sphere radius must be positive.");

  sp_func_t* s = sphere_sp_func_new(x, r, INWARD_NORMAL);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static docstring_t* sphere_doc()
{
  return docstring_from_string("sphere(x, R) - Creates a signed distance function for a sphere centered\n"
                               "  at the point x with the radius R.");
}

static int rect_prism(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_isboundingbox(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "F = rect_prism(bbox)");
  }

  // Get the arguments.
  bbox_t* bbox = lua_toboundingbox(lua, 1);
  sp_func_t* prism = rect_prism_sp_func_from_bbox(bbox);
  lua_pushscalarfunction(lua, st_func_from_sp_func(prism));
  prism = NULL;
  return 1;
}

static docstring_t* rect_prism_doc()
{
  return docstring_from_string("rect_prism(bbox) - Creates a signed distance function for a rectangular\n"
                               "  prism occupying the space within the bounding box bbox.");
}

void interpreter_register_sd_functions(interpreter_t* interp);
void interpreter_register_sd_functions(interpreter_t* interp)
{
  interpreter_register_global_table(interp, "signed_distance_functions", signed_distance_functions_doc());
  interpreter_register_global_method(interp, "signed_distance_functions", "sphere", sphere, sphere_doc());
  interpreter_register_global_method(interp, "signed_distance_functions", "rect_prism", rect_prism, rect_prism_doc());
}

