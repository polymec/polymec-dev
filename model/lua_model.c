// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/lua_model.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static lua_module_function model_funcs[] = {
  {NULL, NULL}
};

static int m_tostring(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_pushfstring(L, "model '%s'", model_name(m));
  return 1;
}

static int m_load(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Argument must be a step to load.");
  int step = (int)(lua_tointeger(L, 2));
  model_load(m, step);
  return 0;
}

static int m_save(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  model_save(m);
  return 0;
}

static int m_advance(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument must be a maximum time step size.");
  real_t dt_max = (real_t)(lua_tonumber(L, 2));
  real_t dt_actual = model_advance(m, dt_max);
  lua_pushnumber(L, dt_actual);
  return 1;
}


static int m_run(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_istable(L, 2))
    luaL_error(L, "Argument must be a table with entries t1, t2, and max_steps.");
  lua_getfield(L, 2, "t1");
  if (!lua_isnumber(L, -1))
    luaL_error(L, "t1 must be a number.");
  real_t t1 = (real_t)(lua_tonumber(L, -1));
  lua_getfield(L, 2, "t2");
  if (!lua_isnumber(L, -1))
    luaL_error(L, "t2 must be a number.");
  real_t t2 = (real_t)(lua_tonumber(L, -1));
  if (t2 < t1)
    luaL_error(L, "t2 must be greater than or equal to t1.");
  lua_getfield(L, 2, "max_steps");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "max_steps must be an integer.");
  int max_steps = (int)(lua_tointeger(L, -1));
  model_run(m, t1, t2, max_steps);

  return 0; // No return value.
}

static lua_class_method model_methods[] = {
  {"load", m_load},
  {"save", m_save},
  {"advance", m_advance},
  {"run", m_run},
  {"__tostring", m_tostring},
  {NULL, NULL}
};

int lua_register_model_modules(lua_State* L)
{
  lua_register_class(L, "model", model_funcs, model_methods);
  return 0;
}

void lua_push_model(lua_State* L, model_t* m)
{
  lua_push_object(L, "model", m, DTOR(model_free));
}

bool lua_is_model(lua_State* L, int index)
{
  return lua_is_object(L, index, "model");
}

model_t* lua_to_model(lua_State* L, int index)
{
  return (model_t*)lua_to_object(L, index, "model");
}

void lua_push_model_probe(lua_State* L, model_probe_t* p)
{
  lua_push_object(L, "model_probe", p, DTOR(model_probe_free));
}

bool lua_is_model_probe(lua_State* L, int index)
{
  return lua_is_object(L, index, "model_probe");
}

model_probe_t* lua_to_model_probe(lua_State* L, int index)
{
  return (model_probe_t*)lua_to_object(L, index, "model_probe");
}

