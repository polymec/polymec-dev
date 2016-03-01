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

static const char* ball_and_jack_usage = 
  "U = scalar_functions.ball_and_jack(x0, ball_radius, ball_value, jack_size, jack_value)\n"
  "  Returns a scalar function that creates a 3D pattern of a jack\n"
  "  of given extent inside a ball of a given radius, centered at the\n"
  "  point x0. Inside the jack, the function adopts the value specified\n"
  "  by jack_value; outside the jack but within the ball, the function\n"
  "  adopts ball_value; outsite that, the function is zero.";

typedef struct
{
  point_t x0;
  real_t ball_radius, ball_value;
  real_t jack_length, jack_thickness, jack_value;
} bnj_t;

static void bnj_eval(void* context, point_t* x, real_t t, real_t* val)
{
  bnj_t* bnj = context;
  vector_t y;
  point_displacement(&bnj->x0, x, &y);
  if (vector_mag(&y) < bnj->ball_radius) // inside the ball
  {
    if (((ABS(y.x) <= 0.5*bnj->jack_length) && 
         (ABS(y.y) <= 0.5*bnj->jack_thickness) && 
         (ABS(y.z) <= 0.5*bnj->jack_thickness)) || 
        ((ABS(y.x) <= 0.5*bnj->jack_thickness) && 
         (ABS(y.y) <= 0.5*bnj->jack_length) && 
         (ABS(y.z) <= 0.5*bnj->jack_thickness)) || 
        ((ABS(y.x) <= 0.5*bnj->jack_thickness) && 
         (ABS(y.y) <= 0.5*bnj->jack_thickness) && 
         (ABS(y.z) <= 0.5*bnj->jack_length)))
      *val = bnj->jack_value;
    else
      *val = bnj->ball_value;
  }
  else
    *val = 0.0;
}

static int ball_and_jack(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 5) || !lua_ispoint(lua, 1) || 
      !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3) ||
      !lua_isnumber(lua, 4) || !lua_isnumber(lua, 5))
    return luaL_error(lua, ball_and_jack_usage);

  // Get the arguments.
  point_t* x0 = lua_topoint(lua, 1);
  real_t ball_radius = (real_t)lua_tonumber(lua, 2);
  if (ball_radius <= 0.0)
    return luaL_error(lua, "ball_and_jack: ball radius must be positive.");
  real_t ball_value = (real_t)lua_tonumber(lua, 3);
  real_t jack_size = (real_t)lua_tonumber(lua, 4);
  real_t jack_thickness = jack_size / 4.0;
  if (jack_size <= 0.0)
    return luaL_error(lua, "ball_and_jack: jack size must be positive.");
  else if (sqrt(0.25*jack_size*jack_size + 0.25*jack_thickness*jack_thickness) > ball_radius)
    return luaL_error(lua, "ball_and_jack: jack size is too large to fit within ball radius.");

  real_t jack_value = (real_t)lua_tonumber(lua, 5);
  bnj_t* bnj = polymec_malloc(sizeof(bnj_t));
  bnj->x0 = *x0;
  bnj->ball_radius = ball_radius;
  bnj->ball_value = ball_value;
  bnj->jack_length = jack_size;
  bnj->jack_thickness = jack_thickness;
  bnj->jack_value = jack_value;
  st_func_vtable vtable = {.eval = bnj_eval, .dtor = polymec_free};
  char bnj_name[1025];
  snprintf(bnj_name, 1024, 
    "ball_and_jack(x0 = (%g, %g, %g), ball_radius = %g, ball_value = %g, jack_size = %g, jack_value = %g)",
    x0->x, x0->y, x0->z, ball_radius, ball_value, jack_size, jack_value);
  st_func_t* func = st_func_new(bnj_name, bnj, vtable, 
                                ST_FUNC_HETEROGENEOUS, ST_FUNC_CONSTANT, 1);
  lua_pushscalarfunction(lua, func);
  return 1;
}

void interpreter_register_scalar_functions(interpreter_t* interp)
{
  interpreter_register_global_table(interp, "scalar_functions", NULL);
  interpreter_register_global_method(interp, "scalar_functions", "ball_and_jack", ball_and_jack, docstring_from_string(ball_and_jack_usage));
}

