// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static const char* box_usage = 
  "I = indicators.box(bounding_box)\n"
  "  Returns an indicator function that returns 1 for points inside the given\n"
  "  bounding box, 0 outside.";

static void box_eval(void* context, point_t* x, real_t t, real_t* val)
{
  bbox_t* B = context;
  if (bbox_contains(B, x))
    *val = 1.0;
  else
    *val = 0.0;
}

static int box(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_isboundingbox(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "F = indicatorsbox(bbox)");
  }

  // Get the arguments.
  bbox_t* bbox = lua_toboundingbox(lua, 1);
  sp_func_t* prism = rect_prism_sp_func_from_bbox(bbox);
  lua_pushscalarfunction(lua, st_func_from_sp_func(prism));
  prism = NULL;
  return 1;
}

static const char* sphere_usage = 
  "I = indicators.sphere(x0, R)\n"
  "  Returns an indicator function that returns 1 for points inside a sphere\n"
  "  of radius R centered at the point x0 and 0 outside that radius.";

typedef struct
{
  point_t x0;
  real_t R;
} sphere_t;

static void sphere_eval(void* context, point_t* x, real_t t, real_t* val)
{
  sphere_t* S = context;
  if (point_distance(x, &S->x0) <= S->R0)
    *val = 1.0;
  else
    *val = 0.0;
}

static int sphere(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispoint(lua, 1) || !lua_isnumber(lua, 2))
  {
    return luaL_error(lua, sphere_usage);
  }

  // Get the arguments.
  point_t* x = lua_topoint(lua, 1);
  real_t r = (real_t)lua_tonumber(lua, 2);
  if (r <= 0.0)
    return luaL_error(lua, "Sphere radius must be positive.");

//  sp_func_t* s = sp_func_new(name, INWARD_NORMAL);
//  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

void interpreter_register_indicators(interpreter_t* interp)
{
  interpreter_register_global_table(interp, "indicators", NULL);
  interpreter_register_global_method(interp, "indicators", "box", box, NULL);
  interpreter_register_global_method(interp, "indicators", "sphere", sphere, NULL);
  interpreter_register_global_method(interp, "indicators", "above_plane", above_plane, NULL);
  interpreter_register_global_method(interp, "indicators", "above_or_on_plane", above_or_on_plane, NULL);
  interpreter_register_global_method(interp, "indicators", "below_plane", below_plane, NULL);
  interpreter_register_global_method(interp, "indicators", "below_or_on_plane", below_or_on_plane, NULL);
}

