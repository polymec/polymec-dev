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

static const char* rigid_body_rotation_usage = 
  "V = vector_functions.rigid_body_rotation{origin = (x0, y0, z0), angular_velocity = (wx, wy, wz)}\n"
  "  Creates a velocity field representing a rigid body rotating in 3D space.\n"
  "  Arguments are:\n"
  "    origin           - Point about which the body rotates.\n"
  "    angular_velocity - Angular velocity (pseudo)vector describing the rotation.";

typedef struct 
{
  point_t x0; // center of rotation
  vector_t omega; // Angular velocity pseudovector.
} rbr_t;

static void rbr_eval(void* context, point_t* x, real_t t, real_t* val)
{
  rbr_t* rbr = context;

  // Displacement vector r.
  vector_t r;
  point_displacement(&rbr->x0, x, &r);

  // V = dr/dt = omega x r.
  vector_cross(&rbr->omega, &r, (vector_t*)val);
}

static int rigid_body_rotation(lua_State* lua)
{
  // Check the number of arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_istable(lua, 1))
    return luaL_error(lua, rigid_body_rotation_usage);

  // Get the origin if it's given.
  lua_pushstring(lua, "origin"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (!lua_isnil(lua, -1) && !lua_ispoint(lua, -1))
    return luaL_error(lua, rigid_body_rotation_usage);
  point_t* x0 = lua_topoint(lua, -1);
  lua_pop(lua, 1);

  // Get the angular velocity vector.
  lua_pushstring(lua, "angular_velocity"); // pushes key onto stack
  lua_gettable(lua, 1); // replaces key with value
  if (lua_isnil(lua, -1) || !lua_isvector(lua, -1))
    return luaL_error(lua, rigid_body_rotation_usage);
  vector_t* omega = lua_tovector(lua, -1);
  lua_pop(lua, 1);

  // Set up the rigid body rotation thingy.
  rbr_t* rbr = polymec_malloc(sizeof(rbr_t));
  if (x0 != NULL)
    rbr->x0 = *x0;
  else
    rbr->x0.x = rbr->x0.y = rbr->x0.z = 0.0;
  rbr->omega = *omega;

  // Create and return the function.
  char name[1025];
  snprintf(name, 1024, "Rigid body rotation(x0 = (%g, %g, %g), omega = (%g, %g, %g))",
           rbr->x0.x, rbr->x0.y, rbr->x0.z, omega->x, omega->y, omega->z);
  st_func_vtable vtable = {.eval = rbr_eval, .dtor = polymec_free};
  st_func_t* func = st_func_new(name, rbr, vtable, ST_FUNC_HETEROGENEOUS, ST_FUNC_CONSTANT, 3);

  lua_pushvectorfunction(lua, func);
  return 1;
}

void interpreter_register_vector_functions(interpreter_t* interp)
{
  interpreter_register_global_table(interp, "vector_functions", NULL);
  interpreter_register_global_method(interp, "vector_functions", "rigid_body_rotation", rigid_body_rotation, docstring_from_string(rigid_body_rotation_usage));
}

