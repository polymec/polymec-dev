// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"
#include "model/lua_model.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

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
    return luaL_error(L, "Argument must be a step to load.");
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
    return luaL_error(L, "Argument must be a maximum time step size.");
  real_t dt_max = lua_to_real(L, 2);
  real_t dt_actual = model_advance(m, dt_max);
  lua_pushnumber(L, dt_actual);
  return 1;
}

static int m_run(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_istable(L, 2))
    return luaL_error(L, "Argument must be a table with entries t1, t2, and max_steps.");
  lua_getfield(L, 2, "t1");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "t1 must be a number.");
  real_t t1 = lua_to_real(L, -1);
  lua_getfield(L, 2, "t2");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "t2 must be a number.");
  real_t t2 = lua_to_real(L, -1);
  if (t2 < t1)
    return luaL_error(L, "t2 must be greater than or equal to t1.");
  lua_getfield(L, 2, "max_steps");
  if (!lua_isinteger(L, -1))
    return luaL_error(L, "max_steps must be an integer.");
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

typedef struct 
{
  lua_State* L;
} lua_probe_t;

static void lp_acquire(void* context, real_t t, int rank, size_t* shape, real_t* data)
{
  lua_probe_t* p = context;

  // Get the callable thingy.
  lua_pushlightuserdata(p->L, p);
  lua_gettable(p->L, LUA_REGISTRYINDEX);

  // If it's an object, call thing thing with itself as the first argument.
  if (lua_istable(p->L, -1))
  {
    lua_pushnumber(p->L, (double)t);
    lua_push_ndarray(p->L, rank, shape, data, LUA_ARRAY_REAL);
    lua_call(p->L, 3, 1);
  }
  // Otherwise, just all it with (t, val).
  else 
  {
    ASSERT(lua_isfunction(p->L, -1));
    lua_pushnumber(p->L, (double)t);
    lua_push_ndarray(p->L, rank, shape, data, LUA_ARRAY_REAL);
    lua_call(p->L, 2, 1);
  }
}

static void lp_dtor(void* context)
{
  lua_probe_t* p = context;

  // Remove the callable thing from our registry.
  lua_pushlightuserdata(p->L, p);
  lua_pushnil(p->L);
  lua_settable(p->L, LUA_REGISTRYINDEX);

  // Kill. 
  polymec_free(p);
}

static int p_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    return luaL_error(L, "Arguments must be 1) a callable object needed and 2) a shape array for the data.");
  if (!lua_istable(L, 1) && !lua_isfunction(L, 1))
    return luaL_error(L, "Argument 1 must be a probe function or a callable object.");
  if (!lua_istable(L, 2))
    return luaL_error(L, "Argument 2 must be a shape array for acquired data.");

  // Jot down the Lua state, and stash the given object in Lua's registry.
  lua_probe_t* p = polymec_malloc(sizeof(lua_probe_t));
  p->L = L;
  lua_pushlightuserdata(L, (void*)p); // key   <-- p
  lua_pushvalue(L, 1);                // value <-- callable object.
  lua_settable(L, LUA_REGISTRYINDEX);

  // Get the shape of the acquired data array.
  int rank = (int)lua_rawlen(L, 2);
  size_t shape[rank];
  for (int i = 1; i <= rank; ++i)
  {
    lua_rawgeti(L, 2, (lua_Integer)i);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Shape array must contain only integers.");
    shape[i-1] = (int)lua_tointeger(L, -1);
  }

  // Set up our probe.
  model_probe_vtable vtable = {.acquire = lp_acquire, .dtor = lp_dtor};
  model_probe_t* probe = model_probe_new("Lua probe", rank, shape, p, vtable);
  lua_push_model_probe(L, probe);
  return 1;
}

static lua_module_function model_probe_funcs[] = {
  {"new", p_new},
  {NULL, NULL}
};

static int p_acquire(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_isnumber(L, 2))
    return luaL_error(L, "Argument must be a time.");

  // Do the acquisition.
  model_probe_t* p = lua_to_model_probe(L, 1);
  ASSERT(p != NULL);
  real_t t = lua_tonumber(L, 2);
  real_t* val = model_probe_new_array(p);
  model_probe_acquire(p, t, val);

  // Return the acquired array as an opaque type.
  lua_pushlightuserdata(L, val);
  return 1;
}

static int p_tostring(lua_State* L)
{
  model_probe_t* p = lua_to_model_probe(L, 1);
  lua_pushfstring(L, "model_probe '%s'", model_probe_name(p));
  return 1;
}

static lua_class_method model_probe_methods[] = {
  {"acquire", p_acquire},
  {"__tostring", p_tostring},
  {NULL, NULL}
};

int lua_register_model_modules(lua_State* L)
{
  lua_register_class(L, "model", NULL, model_methods);
  lua_register_class(L, "model_probe", model_probe_funcs, model_probe_methods);

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

