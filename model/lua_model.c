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
  real_t dt_max = (real_t)(lua_tonumber(L, 2));
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
  real_t t1 = (real_t)(lua_tonumber(L, -1));
  lua_getfield(L, 2, "t2");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "t2 must be a number.");
  real_t t2 = (real_t)(lua_tonumber(L, -1));
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

static void lp_acquire(void* context, real_t t, tensor_t* datum)
{
  lua_probe_t* p = context;

  // Get the callable thingy.
  lua_pushlightuserdata(p->L, p);
  lua_gettable(p->L, LUA_REGISTRYINDEX);

  // If it's an object, call thing thing with itself as the first argument.
  if (lua_istable(p->L, -1))
  {
    lua_pushnumber(p->L, (double)t);
    lua_push_tensor(p->L, datum);
    lua_call(p->L, 3, 1);
  }
  // Otherwise, just all it with (t, val).
  else 
  {
    ASSERT(lua_isfunction(p->L, -1));
    lua_pushnumber(p->L, (double)t);
    lua_push_tensor(p->L, datum);
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
    return luaL_error(L, "Arguments must be 1) a callable object needed and 2) a shape array for tensor data.");
  if (!lua_istable(L, 1) && !lua_isfunction(L, 1))
    return luaL_error(L, "Argument 1 must be a probe function or a callable object.");
  if (!lua_istable(L, 2))
    return luaL_error(L, "Argument 2 must be a shape array for acquired tensor data.");

  // Jot down the Lua state, and stash the given object in Lua's registry.
  lua_probe_t* p = polymec_malloc(sizeof(lua_probe_t));
  p->L = L;
  lua_pushlightuserdata(L, (void*)p); // key   <-- p
  lua_pushvalue(L, 1);                // value <-- callable object.
  lua_settable(L, LUA_REGISTRYINDEX);

  // Get the shape of tensor-valued acquired data.
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
  tensor_t* val = model_probe_new_datum(p);
  model_probe_acquire(p, t, val);

  // Return the acquired tensor.
  lua_push_tensor(L, val);
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

static int c_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args > 0)
    return luaL_error(L, "Too many arguments! None needed.");
  model_data_channel_t* c = local_data_channel_new();
  lua_push_model_data_channel(L, c);
  return 1;
}

static lua_module_function model_data_channel_funcs[] = {
  {"new", c_new},
  {NULL, NULL}
};

static int c_add_output(lua_State* L)
{
  model_data_channel_t* c = lua_to_model_data_channel(L, 1);
  int num_args = lua_gettop(L);
  if (num_args != 3)
    return luaL_error(L, "Arguments must be 1) a list of data names and 2) an output.");

  // This fails if we're not using a local data channel!
  if (strcmp(model_data_channel_name(c), "Local data channel") != 0)
    return luaL_error(L, "Object must be a local_model_data_channel.");

  // We need a list of names.
  if (!lua_istable(L, 2))
    return luaL_error(L, "Argument 1 must be a list of names.");
  string_array_t* data_names = string_array_new();
  size_t len = lua_rawlen(L, 2);
  for (size_t i = 1; i <= len; ++i)
  {
    lua_rawgeti(L, 2, (lua_Integer)i);
    if (!lua_isstring(L, -1))
    {
      string_array_free(data_names);
      return luaL_error(L, "Argument 1 must be a list of names.");
    }
    string_array_append_with_dtor(data_names, 
                                  string_dup(lua_tostring(L, -1)),
                                  string_free);
  }

  // We need a model_local_output.
  if (!lua_is_model_local_output(L, 3))
    return luaL_error(L, "Argument 2 must be a model_local_data_output.");

  // Do it.
  model_local_output_t* o = lua_to_model_local_output(L, 3);
  local_data_channel_add_output(c, data_names, o);

  return 0;
}

static int c_put(lua_State* L)
{
  model_data_channel_t* c = lua_to_model_data_channel(L, 1);
  int num_args = lua_gettop(L);
  if (num_args != 4)
    return luaL_error(L, "Arguments must include 1) a time, 2) a datum name, and 3) a tensor.");
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Argument 1 must be a valid time.");
  if (!lua_isstring(L, 3))
    return luaL_error(L, "Argument 2 must be a valid datum name.");
  if (!lua_is_tensor(L, 4))
    return luaL_error(L, "Argument 3 must be a tensor.");
  real_t t = (real_t)(lua_tonumber(L, 2));
  char* name = (char*)lua_tostring(L, 3);
  tensor_t* val = lua_to_tensor(L, 4);
  model_data_channel_put(c, t, name, val);
  return 0;
}

static int c_tostring(lua_State* L)
{
  model_data_channel_t* c = lua_to_model_data_channel(L, 1);
  lua_pushfstring(L, "model_data_channel '%s'", model_data_channel_name(c));
  return 1;
}

static lua_class_method model_data_channel_methods[] = {
  {"add_output", c_add_output},
  {"put", c_put},
  {"__tostring", c_tostring},
  {NULL, NULL}
};

typedef struct 
{
  lua_State* L;
} lua_output_t;

static void lo_put(void* context, real_t t, char* name, tensor_t* datum)
{
  lua_output_t* o = context;

  // Get the callable thingy.
  lua_pushlightuserdata(o->L, o);
  lua_gettable(o->L, LUA_REGISTRYINDEX);

  // If it's an object, call thing thing with itself as the first argument.
  if (lua_istable(o->L, -1))
  {
    lua_pushnumber(o->L, (double)t);
    lua_pushstring(o->L, (const char*)name);
    lua_push_tensor(o->L, datum);
    lua_call(o->L, 4, 1);
  }
  // Otherwise, just all it with (t, name, val).
  else 
  {
    ASSERT(lua_isfunction(o->L, -1));
    lua_pushnumber(o->L, (double)t);
    lua_pushstring(o->L, (const char*)name);
    lua_push_tensor(o->L, datum);
    lua_call(o->L, 3, 1);
  }
}

static void lo_dtor(void* context)
{
  lua_output_t* o = context;

  // Remove the callable thing from our registry.
  lua_pushlightuserdata(o->L, o);
  lua_pushnil(o->L);
  lua_settable(o->L, LUA_REGISTRYINDEX);

  // Kill. 
  polymec_free(o);
}

static int o_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || (!lua_isfunction(L, 1) && !lua_istable(L, 1)))
    return luaL_error(L, "Argument should be a callable object.");

  // Jot down the Lua state, and stash the given object in Lua's registry.
  lua_output_t* o = polymec_malloc(sizeof(lua_output_t));
  o->L = L;
  lua_pushlightuserdata(L, (void*)o); // key   <-- o
  lua_pushvalue(L, 1);                // value <-- callable object.
  lua_settable(L, LUA_REGISTRYINDEX);

  // Set up our output.
  model_local_output_vtable vtable = {.put = lo_put, .dtor = lo_dtor};
  model_local_output_t* output = model_local_output_new("Lua output", o, vtable);
  lua_push_model_local_output(L, output);
  return 1;
}

static int o_text(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_isstring(L, 1) || !lua_isstring(L, 2))
    return luaL_error(L, "Arguments are a directory and a prefix.");
  const char* dir = lua_tostring(L, 1);
  const char* prefix = lua_tostring(L, 2);
  model_local_output_t* o = text_model_local_output_new(dir, prefix);
  lua_push_model_local_output(L, o);
  return 1;
}

static lua_module_function model_local_output_funcs[] = {
  {"new", o_new},
  {"text", o_text},
  {NULL, NULL}
};

static int o_put(lua_State* L)
{
  model_local_output_t* o = lua_to_model_local_output(L, 1);
  int num_args = lua_gettop(L);
  if (num_args != 4)
    return luaL_error(L, "Arguments must include 1) a time, 2) a datum name, and 3) a tensor.");
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Argument 1 must be a valid time.");
  if (!lua_isstring(L, 3))
    return luaL_error(L, "Argument 2 must be a valid datum name.");
  if (!lua_is_tensor(L, 4))
    return luaL_error(L, "Argument 3 must be a tensor.");
  real_t t = (real_t)(lua_tonumber(L, 2));
  char* name = (char*)lua_tostring(L, 3);
  tensor_t* val = lua_to_tensor(L, 4);
  model_local_output_put(o, t, name, val);
  return 0;
}

static int o_tostring(lua_State* L)
{
  model_local_output_t* o = lua_to_model_local_output(L, 1);
  lua_pushfstring(L, "model_local_output '%s'", model_local_output_name(o));
  return 1;
}

static lua_class_method model_local_output_methods[] = {
  {"put", o_put},
  {"__tostring", o_tostring},
  {NULL, NULL}
};

extern int lua_register_gfx(lua_State* L);

int lua_register_model_modules(lua_State* L)
{
  lua_register_class(L, "model", NULL, model_methods);
  lua_register_class(L, "model_probe", model_probe_funcs, model_probe_methods);
  lua_register_class(L, "model_data_channel", model_data_channel_funcs, model_data_channel_methods);
  lua_register_class(L, "model_local_output", model_local_output_funcs, model_local_output_methods);

  lua_register_gfx(L);

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

void lua_push_model_data_channel(lua_State* L, model_data_channel_t* c)
{
  lua_push_object(L, "model_data_channel", c, DTOR(model_data_channel_free));
}

bool lua_is_model_data_channel(lua_State* L, int index)
{
  return lua_is_object(L, index, "model_data_channel");
}

model_data_channel_t* lua_to_model_data_channel(lua_State* L, int index)
{
  return (model_data_channel_t*)lua_to_object(L, index, "model_data_channel");
}

void lua_push_model_local_output(lua_State* L, model_local_output_t* o)
{
  lua_push_object(L, "model_local_output", o, DTOR(model_local_output_free));
}

bool lua_is_model_local_output(lua_State* L, int index)
{
  return lua_is_object(L, index, "model_local_output");
}

model_local_output_t* lua_to_model_local_output(lua_State* L, int index)
{
  return (model_local_output_t*)lua_to_object(L, index, "model_local_output");
}

