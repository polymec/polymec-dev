// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/plane_sp_func.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static const char* box_usage = 
  "I = indicators.box(bounding_box)\n"
  "  Returns an indicator function that returns 1 for points inside the given\n"
  "  bounding box, 0 outside.";

static void box_eval(void* context, point_t* x, real_t* val)
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
    return luaL_error(lua, box_usage);

  // Get the arguments.
  bbox_t* B = lua_toboundingbox(lua, 1);
  sp_func_vtable vtable = {.eval = box_eval};
  char name[1025];
  snprintf(name, 1024, "Box indicator (x1 = %g, x2 = %g, y1 = %g, y2 = %g, z1 = %g, z2 = %g)",
           B->x1, B->x2, B->y1, B->y2, B->z1, B->z2);
  sp_func_t* s = sp_func_new(name, B, vtable, SP_FUNC_HETEROGENEOUS, 1);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
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

static void sphere_eval(void* context, point_t* x, real_t* val)
{
  sphere_t* S = context;
  if (point_distance(x, &S->x0) <= S->R)
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
  real_t R = (real_t)lua_tonumber(lua, 2);
  if (R <= 0.0)
    return luaL_error(lua, "Sphere radius must be positive.");

  sphere_t* S = polymec_malloc(sizeof(sphere_t));
  S->x0 = *x;
  S->R = R;
  sp_func_vtable vtable = {.eval = sphere_eval, .dtor = polymec_free};
  char name[1025];
  snprintf(name, 1024, "Sphere indicator (x0 = (%g, %g, %g), R = %g)", x->x, x->y, x->z, R);
  sp_func_t* s = sp_func_new(name, S, vtable, SP_FUNC_HETEROGENEOUS, 1);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static const char* above_plane_usage = 
  "I = indicators.above_plane(n, x)\n"
  "  Returns an indicator function that returns 1 for points above the\n"
  "  plane with normal vector n containing point x, and 0 in or below the plane.\n"
  "  A point is 'above' the plane if the normal vector points toward the half\n"
  "  space containing the point.";

static void above_plane_eval(void* context, point_t* x, real_t* val)
{
  sp_func_t* P = context;
  real_t D;
  sp_func_eval(P, x, &D);
  if (D < 0.0)
    *val = 1.0;
  else
    *val = 0.0;
}

static int above_plane(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_isvector(lua, 1) || !lua_ispoint(lua, 2))
  {
    return luaL_error(lua, above_plane_usage);
  }

  // Get the arguments.
  vector_t* n = lua_tovector(lua, 1);
  point_t* x = lua_topoint(lua, 2);

  sp_func_t* P = plane_sp_func_new(n, x);
  sp_func_vtable vtable = {.eval = above_plane_eval};
  char name[1025];
  snprintf(name, 1024, "Above-plane indicator (x = (%g, %g, %g), n = (%g, %g, %g))", 
           x->x, x->y, x->z, n->x, n->y, n->z);
  sp_func_t* s = sp_func_new(name, P, vtable, SP_FUNC_HETEROGENEOUS, 1);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static const char* above_or_in_plane_usage = 
  "I = indicators.above_or_in_plane(n, x)\n"
  "  Returns an indicator function that returns 1 for points above or in the\n"
  "  plane with normal vector n containing point x, and 0 below the plane.\n"
  "  A point is 'above' the plane if the normal vector points toward the half\n"
  "  space containing the point.";

static void above_or_in_plane_eval(void* context, point_t* x, real_t* val)
{
  sp_func_t* P = context;
  real_t D;
  sp_func_eval(P, x, &D);
  if (D <= 0.0)
    *val = 1.0;
  else
    *val = 0.0;
}

static int above_or_in_plane(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_isvector(lua, 1) || !lua_ispoint(lua, 2))
  {
    return luaL_error(lua, above_or_in_plane_usage);
  }

  // Get the arguments.
  vector_t* n = lua_tovector(lua, 1);
  point_t* x = lua_topoint(lua, 2);

  sp_func_t* P = plane_sp_func_new(n, x);
  sp_func_vtable vtable = {.eval = above_or_in_plane_eval};
  char name[1025];
  snprintf(name, 1024, "Above-or-in-plane indicator (x = (%g, %g, %g), n = (%g, %g, %g))", 
           x->x, x->y, x->z, n->x, n->y, n->z);
  sp_func_t* s = sp_func_new(name, P, vtable, SP_FUNC_HETEROGENEOUS, 1);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static const char* below_plane_usage = 
  "I = indicators.below_plane(n, x)\n"
  "  Returns an indicator function that returns 1 for points below the\n"
  "  plane with normal vector n containing point x, and 0 in or above the plane.\n"
  "  A point is 'below' the plane if the normal vector points away from the half\n"
  "  space containing the point.";

static void below_plane_eval(void* context, point_t* x, real_t* val)
{
  sp_func_t* P = context;
  real_t D;
  sp_func_eval(P, x, &D);
  if (D > 0.0)
    *val = 1.0;
  else
    *val = 0.0;
}

static int below_plane(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_isvector(lua, 1) || !lua_ispoint(lua, 2))
  {
    return luaL_error(lua, below_plane_usage);
  }

  // Get the arguments.
  vector_t* n = lua_tovector(lua, 1);
  point_t* x = lua_topoint(lua, 2);

  sp_func_t* P = plane_sp_func_new(n, x);
  sp_func_vtable vtable = {.eval = below_plane_eval};
  char name[1025];
  snprintf(name, 1024, "Below-plane indicator (x = (%g, %g, %g), n = (%g, %g, %g))", 
           x->x, x->y, x->z, n->x, n->y, n->z);
  sp_func_t* s = sp_func_new(name, P, vtable, SP_FUNC_HETEROGENEOUS, 1);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static const char* below_or_in_plane_usage = 
  "I = indicators.below_or_in_plane(n, x)\n"
  "  Returns an indicator function that returns 1 for points below or in the\n"
  "  plane with normal vector n containing point x, and 0 above the plane.\n"
  "  A point is 'below' the plane if the normal vector points away from the half\n"
  "  space containing the point.";

static void below_or_in_plane_eval(void* context, point_t* x, real_t* val)
{
  sp_func_t* P = context;
  real_t D;
  sp_func_eval(P, x, &D);
  if (D <= 0.0)
    *val = 1.0;
  else
    *val = 0.0;
}

static int below_or_in_plane(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_isvector(lua, 1) || !lua_ispoint(lua, 2))
  {
    return luaL_error(lua, below_or_in_plane_usage);
  }

  // Get the arguments.
  vector_t* n = lua_tovector(lua, 1);
  point_t* x = lua_topoint(lua, 2);

  sp_func_t* P = plane_sp_func_new(n, x);
  sp_func_vtable vtable = {.eval = below_or_in_plane_eval};
  char name[1025];
  snprintf(name, 1024, "Below-or-in-plane indicator (x = (%g, %g, %g), n = (%g, %g, %g))", 
           x->x, x->y, x->z, n->x, n->y, n->z);
  sp_func_t* s = sp_func_new(name, P, vtable, SP_FUNC_HETEROGENEOUS, 1);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static const char* indicators_doc = 
  "indicators -- Indicator functions that return 1 inside of a region and 0 outside.";

void interpreter_register_indicators(interpreter_t* interp)
{
  interpreter_register_global_table(interp, "indicators", docstring_from_string(indicators_doc));
  interpreter_register_global_method(interp, "indicators", "box", box, docstring_from_string(box_usage));
  interpreter_register_global_method(interp, "indicators", "sphere", sphere, docstring_from_string(sphere_usage));
  interpreter_register_global_method(interp, "indicators", "above_plane", above_plane, docstring_from_string(above_plane_usage));
  interpreter_register_global_method(interp, "indicators", "above_or_in_plane", above_or_in_plane, docstring_from_string(above_or_in_plane_usage));
  interpreter_register_global_method(interp, "indicators", "below_plane", below_plane, docstring_from_string(below_plane_usage));
  interpreter_register_global_method(interp, "indicators", "below_or_in_plane", below_or_in_plane, docstring_from_string(below_or_in_plane_usage));
}

