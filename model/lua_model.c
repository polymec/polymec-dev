// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"
#include "core/array.h"
#include "core/declare_nd_array.h"
#include "model/lua_model.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

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

static int m_add_probe(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_is_probe(L, 2))
    return luaL_error(L, "Argument 2 must be a probe.");
  if (!lua_istable(L, 3))
    return luaL_error(L, "Argument 3 must be a table of acquisition times.");

  // Build an array of acquisition times.
  real_array_t* times = real_array_new();
  lua_pushnil(L);
  while (lua_next(L, 3))
  {
    // Key is at index -2, value is at -1.
    if (!lua_isinteger(L, -2))
      luaL_error(L, "Argument 3 must be a array of acquisition times.");
    if (!lua_isnumber(L, -1))
      luaL_error(L, "Argument 3 must be a array of acquisition times.");
    real_array_append(times, lua_to_real(L, -1));
    lua_pop(L, 1);
  }

  // Add the probe.
  probe_t* probe = lua_to_probe(L, 2);
  model_add_probe(m, probe, times->data, times->size);
  real_array_release_data_and_free(times);

  return 0;
}

static int m_tostring(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_pushfstring(L, "model '%s'", model_name(m));
  return 1;
}

static lua_class_method model_methods[] = {
  {"load", m_load},
  {"save", m_save},
  {"advance", m_advance},
  {"run", m_run},
  {"add_probe", m_add_probe},
  {"__tostring", m_tostring},
  {NULL, NULL}
};

typedef struct 
{
  lua_State* L;
} lua_probe_t;

static void table_to_probe_data(lua_State* L, int index, probe_data_t* data)
{
  int rank = data->rank;
  ASSERT((rank >= 0) && (rank <= 4));
  if (rank == 0)
  {
    if (!lua_isnumber(L, index))
      luaL_error(L, "probe:acquire: data is not a number.");
    else
      *data->data = lua_to_real(L, index);
  }
  else if (lua_istable(L, index))
  {
    for (size_t i = 1; i <= data->shape[0]; ++i)
    {
      lua_rawgeti(L, index, (int)i);
      if (rank == 1)
      {
        if (lua_isnumber(L, -1))
          data->data[i-1] = lua_to_real(L, -1);
        else if (lua_isnil(L, -1))
          data->data[i-1] = 0.0;
        else
          luaL_error(L, "probe:acquire: data[%d] is not a number.", (int)i);
      }
      else 
      {
        for (size_t j = 1; j <= data->shape[1]; ++j)
        {
          lua_rawgeti(L, -1, (int)j);
          if (rank == 2)
          {
            DECLARE_2D_ARRAY(real_t, xij, data->data, data->shape[0], data->shape[1]);
            if (lua_isnumber(L, -1))
              xij[i-1][j-1] = lua_to_real(L, -1);
            else if (lua_isnil(L, -1))
              xij[i-1][j-1] = 0.0;
            else
              luaL_error(L, "probe:acquire: data[%d][%d] is not a number.", (int)i, (int)j);
          }
          else
          {
            for (size_t k = 1; k <= data->shape[2]; ++k)
            {
              lua_rawgeti(L, -1, (int)k);
              if (rank == 3)
              {
                DECLARE_3D_ARRAY(real_t, xijk, data->data, data->shape[0], data->shape[1], data->shape[2]);
                if (lua_isnumber(L, -1))
                  xijk[i-1][j-1][k-1] = lua_to_real(L, -1);
                else if (lua_isnil(L, -1))
                  xijk[i-1][j-1][k-1] = 0.0;
                else
                  luaL_error(L, "probe:acquire: data[%d][%d][%d] is not a number.", (int)i, (int)j, (int)k);
              }
              else // rank == 4
              {
                for (size_t l = 1; l <= data->shape[3]; ++l)
                {
                  lua_rawgeti(L, -1, (int)l);
                  DECLARE_4D_ARRAY(real_t, xijkl, data->data, data->shape[0], data->shape[1], data->shape[2], data->shape[3]);
                  if (lua_isnumber(L, -1))
                    xijkl[i-1][j-1][k-1][l-1] = lua_to_real(L, -1);
                  else if (lua_isnil(L, -1))
                    xijkl[i-1][j-1][k-1][l-1] = 0.0;
                  else
                    luaL_error(L, "probe:acquire: data[%d][%d][%d] is not a number.", (int)i, (int)j, (int)k);
                  lua_pop(L, 1); // pop l index
                }
              }
              lua_pop(L, 1); // pop k index
            }
          }
          lua_pop(L, 1); // pop j index
        }
      }
      lua_pop(L, 1); // pop i index
    }
  }
  else
    luaL_error(L, "probe:acquire: data is not a table.");
}

static void lp_acquire(void* context, real_t t, probe_data_t* data)
{
  lua_probe_t* p = context;

  // Get the callable thingy.
  lua_pushlightuserdata(p->L, p);
  lua_gettable(p->L, LUA_REGISTRYINDEX);

  // If it's an object, call thing thing with itself as the first argument.
  if (lua_istable(p->L, -1))
  {
    lua_push_real(p->L, t);
    lua_call(p->L, 2, 1);
  }
  // Otherwise, just call it with t.
  else 
  {
    ASSERT(lua_isfunction(p->L, -1));
    lua_push_real(p->L, t);
    lua_call(p->L, 1, 1);
  }

  // The return value is a table, which we use to fill our probe_data.
  table_to_probe_data(p->L, -1, data);
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
  if (num_args != 3)
    return luaL_error(L, "Arguments must be 1) a quantity name, 2) a callable object, and 3) a shape array for the data.");
  if (!lua_isstring(L, 1))
    return luaL_error(L, "Argument 1 must be a string.");
  if (!lua_istable(L, 2) && !lua_isfunction(L, 2))
    return luaL_error(L, "Argument 2 must be a probe function or a callable object.");
  if (!lua_istable(L, 3))
    return luaL_error(L, "Argument 3 must be a shape array for acquired data.");

  const char* data_name = lua_tostring(L, 1);

  // Jot down the Lua state, and stash the given object in Lua's registry.
  lua_probe_t* p = polymec_malloc(sizeof(lua_probe_t));
  p->L = L;
  lua_pushlightuserdata(L, (void*)p); // key   <-- p
  lua_pushvalue(L, 2);                // value <-- callable object.
  lua_settable(L, LUA_REGISTRYINDEX);

  // Get the shape of the acquired data array.
  int rank = (int)lua_rawlen(L, 3);
  if (rank > 4)
    return luaL_error(L, "Shape array must be rank 4 or less.");
  size_t shape[rank];
  for (int i = 1; i <= rank; ++i)
  {
    lua_rawgeti(L, 3, (lua_Integer)i);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Shape array must contain only integers.");
    shape[i-1] = (size_t)lua_tointeger(L, -1);
  }

  // Set up our probe.
  probe_vtable vtable = {.acquire = lp_acquire, .dtor = lp_dtor};
  probe_t* probe = probe_new("Lua probe", data_name, rank, shape, p, vtable);
  lua_push_probe(L, probe);
  return 1;
}

static lua_module_function probe_funcs[] = {
  {"new", p_new},
  {NULL, NULL}
};

static int lua_push_probe_data(lua_State* L, probe_data_t* data)
{
  int rank = data->rank;
  if (rank > 4)
    luaL_error(L, "probe:acquire: rank must be 4 or less.");

  if (rank == 0)
    lua_push_real(L, data->data[0]);
  else
  {
    lua_newtable(L);
    for (size_t i = 1; i <= data->shape[0]; ++i)
    {
      if (rank == 1)
      {
        lua_push_real(L, data->data[i-1]);
        lua_rawseti(L, -2, (int)i);
      }
      else
      {
        lua_newtable(L);
        for (size_t j = 1; j <= data->shape[1]; ++j)
        {
          if (rank == 2)
          {
            DECLARE_2D_ARRAY(real_t, xij, data->data, data->shape[0], data->shape[1]);
            lua_push_real(L, xij[i-1][j-1]);
            lua_rawseti(L, -2, (int)j);
          }
          else
          {
            lua_newtable(L);
            for (size_t k = 1; k <= data->shape[2]; ++k)
            {
              if (rank == 3)
              {
                DECLARE_3D_ARRAY(real_t, xijk, data->data, data->shape[0], data->shape[1], data->shape[2]);
                lua_push_real(L, xijk[i-1][j-1][k-1]);
                lua_rawseti(L, -2, (int)k);
              }
              else
              {
                lua_newtable(L); 
                for (size_t l = 1; l <= data->shape[3]; ++l)
                {
                  // rank == 4
                  DECLARE_4D_ARRAY(real_t, xijkl, data->data, data->shape[0], data->shape[1], data->shape[2], data->shape[3]);
                  lua_push_real(L, xijkl[i-1][j-1][k-1][l-1]);
                  lua_rawseti(L, -2, (int)l);
                }
                lua_rawseti(L, -2, (int)k); // push the new table into place
              }
            }
            lua_rawseti(L, -2, (int)j); // push the table into place
          }
        }
        lua_rawseti(L, -2, (int)i); // push the table into place
      }
    }
  }
  return 1;
}

static int p_acquire(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_isnumber(L, 2))
    return luaL_error(L, "Argument must be a time.");

  // Do the acquisition.
  probe_t* p = lua_to_probe(L, 1);
  ASSERT(p != NULL);
  real_t t = lua_to_real(L, 2);
  probe_data_t* data = probe_acquire(p, t);

  // Return the acquired data as a table.
  lua_newtable(L);
  lua_push_real(L, t);
  lua_setfield(L, -2, "time");
  lua_push_probe_data(L, data);
  lua_setfield(L, -2, "data");
  return 1;
}

static int p_tostring(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  lua_pushfstring(L, "probe '%s'", probe_name(p));
  return 1;
}

static lua_class_method probe_methods[] = {
  {"acquire", p_acquire},
  {"__tostring", p_tostring},
  {NULL, NULL}
};

int lua_register_model_modules(lua_State* L)
{
  lua_register_class(L, "model", NULL, model_methods);
  lua_register_class(L, "probe", probe_funcs, probe_methods);

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

void lua_push_probe(lua_State* L, probe_t* p)
{
  lua_push_object(L, "probe", p, DTOR(probe_free));
}

bool lua_is_probe(lua_State* L, int index)
{
  return lua_is_object(L, index, "probe");
}

probe_t* lua_to_probe(lua_State* L, int index)
{
  return (probe_t*)lua_to_object(L, index, "probe");
}

