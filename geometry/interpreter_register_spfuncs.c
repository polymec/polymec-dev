// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/interpreter.h"
#include "geometry/sphere.h"
#include "geometry/rect_prism.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int sphere(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispoint(lua, 1) || !lua_isnumber(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "F = sphere(x, r)");
  }

  // Get the arguments.
  point_t* x = lua_topoint(lua, 1);
  real_t r = (real_t)lua_tonumber(lua, 2);
  if (r <= 0.0)
    return luaL_error(lua, "Sphere radius must be positive.");

  sp_func_t* s = sphere_new(x, r, INWARD_NORMAL);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
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
  sp_func_t* prism = rect_prism_from_bbox(bbox);
  lua_pushscalarfunction(lua, st_func_from_sp_func(prism));
  prism = NULL;
  return 1;
}

void interpreter_register_spfuncs(interpreter_t* interp)
{
  interpreter_register_function(interp, "sphere", sphere, NULL);
  interpreter_register_function(interp, "rect_prism", rect_prism, NULL);
}

