// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

typedef struct
{
  lua_State* L;
} lua_model_t;

static void lm_init(void* context, real_t t)
{
  // Fetch our Lua table from the registry.
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Call the init method.
  lua_getfield(lm->L, -1, "init");
  lua_pushvalue(lm->L, -2);
  lua_pushnumber(lm->L, t);
  lua_call(lm->L, 2, 0);
}

static real_t lm_max_dt(void* context, real_t t, char* reason)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Call the max_dt method.
  lua_getfield(lm->L, -1, "max_dt");
  lua_pushvalue(lm->L, -2);
  lua_pushnumber(lm->L, t);
  lua_call(lm->L, 2, 2);

  if (!lua_isnumber(lm->L, -2))
  {
    luaL_error(lm->L, "max_dt function did not return a dt.");
    return 0.0;
  }
  if (!lua_isstring(lm->L, -1))
  {
    luaL_error(lm->L, "max_dt function did not return a string.");
    return 0.0;
  }
  const char* r = lua_tostring(lm->L, -1);
  strcpy(reason, r);
  return lua_to_real(lm->L, -2);
}

static real_t lm_advance(void* context, real_t max_dt, real_t t)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Do it.
  lua_getfield(lm->L, -1, "advance");
  lua_pushvalue(lm->L, -2);
  lua_pushnumber(lm->L, max_dt);
  lua_pushnumber(lm->L, t);
  lua_call(lm->L, 3, 1);

  if (!lua_isstring(lm->L, -1))
    return luaL_error(lm->L, "advance function did not return a dt.");
  return lua_to_real(lm->L, -1);
}

static bool lm_load(void* context, const char* file_prefix, const char* directory, real_t* time, int step)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Do it!
  lua_getfield(lm->L, -1, "load");
  lua_pushvalue(lm->L, -2);
  lua_pushstring(lm->L, file_prefix);
  lua_pushstring(lm->L, directory);
  lua_pushinteger(lm->L, step);
  lua_call(lm->L, 4, 1);

  if (lua_isnumber(lm->L, -1))
  {
    *time = lua_to_real(lm->L, -1);
    return true;
  }
  return false;
}

static void lm_save(void* context, const char* file_prefix, const char* directory, real_t time, int step)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Call it!
  lua_getfield(lm->L, -1, "save");
  lua_pushvalue(lm->L, -2);
  lua_pushstring(lm->L, file_prefix);
  lua_pushstring(lm->L, directory);
  lua_pushnumber(lm->L, time);
  lua_pushinteger(lm->L, step);
  lua_call(lm->L, 5, 0);
}

static void lm_plot(void* context, const char* file_prefix, const char* directory, real_t time, int step)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Call it!
  lua_getfield(lm->L, -1, "plot");
  lua_pushvalue(lm->L, -2);
  lua_pushstring(lm->L, file_prefix);
  lua_pushstring(lm->L, directory);
  lua_pushnumber(lm->L, time);
  lua_pushinteger(lm->L, step);
  lua_call(lm->L, 5, 0);
}

static void lm_finalize(void* context, int step, real_t t)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);

  // Call it!
  lua_getfield(lm->L, -1, "finalize");
  lua_pushvalue(lm->L, -2);
  lua_pushnumber(lm->L, step);
  lua_pushnumber(lm->L, t);
  lua_call(lm->L, 3, 0);
}

static void lm_dtor(void* context)
{
  lua_model_t* lm = context;
  lua_pushlightuserdata(lm->L, lm);
  lua_gettable(lm->L, LUA_REGISTRYINDEX);
  lua_pushnil(lm->L);
  lua_settable(lm->L, LUA_REGISTRYINDEX);
  polymec_free(context);
}

static int m_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with keys for its name and methods.");

  // Create a table at the top of the stack. This will be the context for
  // the methods.
  lua_newtable(L);
  int obj_index = num_args + 1;

  // Fetch keys and fill in blanks.
  char* name = NULL;
  model_vtable vtable = {.dtor = lm_dtor};
  lua_pushnil(L);
  while (lua_next(L, 1))
  {
    // Key is at index -2, value is at -1.
    if (!lua_isstring(L, -2))
      luaL_error(L, "Keys in table must be strings.");
    const char* key = lua_tostring(L, -2);
    if (strcmp(key, "name") == 0)
    {
      if (!lua_isstring(L, -1))
        luaL_error(L, "Name must be a string.");
      name = string_dup(lua_tostring(L, -1));
      lua_pop(L, 1);
    }
    else
    {
      if (!lua_isfunction(L, -1))
        luaL_error(L, "Methods must be functions.");
      if (strcmp(key, "init") == 0)
      {
        lua_setfield(L, obj_index, "init");
        vtable.init = lm_init;
      }
      else if (strcmp(key, "max_dt") == 0)
      {
        lua_setfield(L, obj_index, "max_dt");
        vtable.max_dt = lm_max_dt;
      }
      else if (strcmp(key, "advance") == 0)
      {
        lua_setfield(L, obj_index, "advance");
        vtable.advance = lm_advance;
      }
      else if (strcmp(key, "load") == 0)
      {
        lua_setfield(L, obj_index, "load");
        vtable.load = lm_load;
      }
      else if (strcmp(key, "save") == 0)
      {
        lua_setfield(L, obj_index, "save");
        vtable.save = lm_save;
      }
      else if (strcmp(key, "plot") == 0)
      {
        lua_setfield(L, obj_index, "plot");
        vtable.plot = lm_plot;
      }
      else if (strcmp(key, "finalize") == 0)
      {
        lua_setfield(L, obj_index, "finalize");
        vtable.finalize = lm_finalize;
      }
    }
  }

  // What are we missing?
  if (vtable.init == NULL)
    luaL_error(L, "init method must be provided.");
  if (vtable.advance == NULL)
    luaL_error(L, "advance method must be provided.");

  // Set up a context to store the Lua state.
  lua_model_t* lm = polymec_malloc(sizeof(lua_model_t));
  lm->L = L;

  // Store the table representing our object in the registry, with
  // context as a key.
  lua_rawsetp(L, LUA_REGISTRYINDEX, lm);

  // Set up the model and get on with it.
  model_t* m = model_new(name, lm, vtable, MODEL_SERIAL);
  lua_push_model(L, m);

  // Clean up.
  string_free(name);

  return 1;
}

static lua_module_function model_funcs[] = {
  {"new", m_new, "model.new{name = NAME, init = INIT, max_dt = MAX_DT, advance = ADVANCE, load = LOAD, save = SAVE, plot = PLOT, finalize = FINALIZE} -> New model with behavior defined by the given functions."},
  {NULL, NULL, NULL}
};

static int m_get_name(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_pushstring(L, model_name(m));
  return 1;
}

static int m_get_context(lua_State* L)
{
  // Fetch our context from the Lua table.
  // It's nil if it's not implemented in Lua.
  model_t* m = lua_to_model(L, 1);
  lua_pushlightuserdata(L, model_context(m));
  lua_gettable(L, LUA_REGISTRYINDEX);
  return 1;
}
static int m_get_step(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_pushinteger(L, model_step(m));
  return 1;
}

static int m_get_time(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_push_real(L, model_time(m));
  return 1;
}

static int m_get_min_dt(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_push_real(L, model_min_dt(m));
  return 1;
}

static int m_set_min_dt(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_is_real(L, 2))
    return luaL_error(L, "Argument must be a positive time interval.");
  real_t dt = lua_to_real(L, 2);
  if (dt <= 0.0)
    return luaL_error(L, "Argument must be a positive time interval.");
  char* reason;
  if (dt > model_max_dt(m, &reason))
    return luaL_error(L, "Model min_dt must not exceed max_dt.");
  model_set_min_dt(m, dt);
  return 0;
}

static int m_get_max_dt(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  char* reason;
  lua_push_real(L, model_max_dt(m, &reason));
  lua_pushstring(L, reason);
  return 2;
}

static int m_set_max_dt(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_is_real(L, 2))
    return luaL_error(L, "Argument must be a positive time interval.");
  real_t dt = lua_to_real(L, 2);
  if (dt <= 0.0)
    return luaL_error(L, "Argument must be a positive time interval.");
  if (dt < model_min_dt(m))
    return luaL_error(L, "Model max_dt must not be smaller than min_dt.");
  model_set_max_dt(m, dt);
  return 0;
}

static int m_get_initial_dt(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  lua_push_real(L, model_initial_dt(m));
  return 1;
}

static int m_set_initial_dt(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (!lua_is_real(L, 2))
    return luaL_error(L, "Argument must be a positive time interval.");
  real_t dt = lua_to_real(L, 2);
  if (dt <= 0.0)
    return luaL_error(L, "Argument must be a positive time interval.");
  model_set_initial_dt(m, dt);
  return 0;
}

static lua_class_field model_fields[] = {
  {"name", m_get_name, NULL},
  {"context", m_get_context, NULL},
  {"step", m_get_step, NULL},
  {"time", m_get_time, NULL},
  {"min_dt", m_get_min_dt, m_set_min_dt},
  {"max_dt", m_get_max_dt, m_set_max_dt},
  {"initial_dt", m_get_initial_dt, m_set_initial_dt},
  {NULL, NULL, NULL}
};

static int m_init(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Argument must be a time at which to initialize the model.");
  real_t t0 = (int)(lua_to_real(L, 2));
  model_init(m, t0);
  return 0;
}

static int m_load(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  if (!lua_isinteger(L, 2))
    return luaL_error(L, "Argument must be a step to load.");
  int step = (int)(lua_tointeger(L, 2));
  bool loaded = model_load(m, step);
  if (!loaded)
    return luaL_error(L, "Model could not be loaded from step %d.", step);
  return 0;
}

static int m_save(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  model_save(m);
  return 0;
}

static int m_plot(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  model_plot(m);
  return 0;
}

static int m_advance(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
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
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  if (!lua_istable(L, 2))
    return luaL_error(L, "Argument must be a table with entries t1, t2, and max_steps.");

  lua_getfield(L, 2, "t2");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "t2 must be a time.");
  real_t t2 = lua_to_real(L, -1);

  real_t t1 = t2;
  lua_getfield(L, 2, "t1");
  if (lua_isnumber(L, -1))
    t1 = lua_to_real(L, -1);
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "t1 must be a time.");

  if (t2 < t1)
    return luaL_error(L, "t2 must be greater than or equal to t1.");

  lua_getfield(L, 2, "max_steps");
  if (!lua_isinteger(L, -1))
    return luaL_error(L, "max_steps must be a number of steps.");
  int max_steps = (int)(lua_tointeger(L, -1));

  lua_getfield(L, 2, "prefix");
  if (lua_isstring(L, -1))
    model_set_prefix(m, lua_tostring(L, -1));
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "prefix must be a string.");

  lua_getfield(L, 2, "directory");
  if (lua_isstring(L, -1))
    model_set_directory(m, lua_tostring(L, -1));
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "directory must be a string.");

  lua_getfield(L, 2, "plot_every");
  if (lua_isnumber(L, -1))
  {
    real_t plot = lua_to_real(L, -1);
    if (plot <= 0.0)
      return luaL_error(L, "plot_every must be positive.");
    model_plot_every(m, plot);
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "plot_every must be a simulation time interval.");

  lua_getfield(L, 2, "save_every");
  if (lua_isinteger(L, -1))
  {
    int save = (int)lua_tointeger(L, -1);
    if (save <= 0)
      return luaL_error(L, "save_every must be positive.");
    model_save_every(m, save);
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "save_every must be a number of steps.");

  int load_step = -1;
  lua_getfield(L, 2, "load_step");
  if (lua_isinteger(L, -1))
  {
    load_step = (int)lua_tointeger(L, -1);
    if (load_step < 0)
      return luaL_error(L, "load_step must be non-negative.");
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "load_step must be a step.");

  model_diag_mode_t diag_mode = MODEL_DIAG_NEAREST_STEP;
  lua_getfield(L, 2, "diag_mode");
  if (lua_isstring(L, -1))
  {
    const char* mode_str = lua_tostring(L, -1);
    if (strcasecmp(mode_str, "exact_time") == 0)
      diag_mode = MODEL_DIAG_EXACT_TIME;
    else if (strcasecmp(mode_str, "nearest_step") != 0)
      return luaL_error(L, "diag_mode must be a 'nearest_step' or 'exact_time'.");
    model_set_diagnostic_mode(m, diag_mode);
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "diag_mode must be a 'nearest_step' or 'exact_time'.");

  real_t min_dt = 0.0;
  lua_getfield(L, 2, "min_dt");
  if (lua_isnumber(L, -1))
  {
    min_dt = lua_to_real(L, -1);
    if (min_dt <= 0.0)
      return luaL_error(L, "min_dt must be positive.");
    model_set_min_dt(m, min_dt);
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "min_dt must be a simulation time interval.");

  real_t max_dt;
  lua_getfield(L, 2, "max_dt");
  if (lua_isnumber(L, -1))
  {
    max_dt = lua_to_real(L, -1);
    if (max_dt < min_dt)
      return luaL_error(L, "max_dt must be positive and greater than or equal to min_dt.");
    model_set_max_dt(m, max_dt);
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "max_dt must be a simulation time interval.");

  // Validate.
  if ((t1 < t2) && (load_step != -1))
    return luaL_error(L, "Only t1 or load_step may be given.");

  if (load_step != -1)
    model_load_from(m, load_step);

  model_run(m, t1, t2, max_steps);

  return 0; // No return value.
}

static int m_handle_signals(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  if (!lua_isboolean(L, 2))
    return luaL_error(L, "Argument 1 must be a boolean flag.");
  bool flag = lua_toboolean(L, 2);
  model_handle_signals(m, flag);
  return 0;
}

static int m_add_probe(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  if (!lua_is_probe(L, 2))
    return luaL_error(L, "Argument 1 must be a probe.");
  if (!lua_istable(L, 3))
    return luaL_error(L, "Argument 2 must be a table of acquisition times.");

  // Build an array of acquisition times.
  real_array_t* times = real_array_new();
  lua_pushnil(L);
  while (lua_next(L, 3))
  {
    // Key is at index -2, value is at -1.
    if (!lua_isinteger(L, -2))
      luaL_error(L, "Argument 2 must be a array of acquisition times.");
    if (!lua_isnumber(L, -1))
      luaL_error(L, "Argument 2 must be a array of acquisition times.");
    real_array_append(times, lua_to_real(L, -1));
    lua_pop(L, 1);
  }

  // Add the probe, transferring its ownership to the model.
  probe_t* probe = lua_to_probe(L, 2);
  lua_transfer_object(L, 2, "probe", LUA_OWNED_BY_C);
  model_add_probe(m, probe, times->data, times->size);
  real_array_free(times);

  return 0;
}

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

static int m_data(lua_State* L)
{
  int num_args = lua_gettop(L);
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  if (num_args == 1) // allll the data.
  {
    int pos = 0;
    char* quantity;
    probe_data_array_t* array;
    lua_newtable(L);
    while (model_next_probe_data(m, &pos, &quantity, &array))
    {
      // We return a table of tables, each with two fields: times and values.
      lua_newtable(L);

      // Times table.
      lua_newtable(L);
      for (size_t i = 1; i <= array->size; ++i)
      {
        lua_push_real(L, array->data[i-1]->time);
        lua_rawseti(L, -2, i);
      }
      lua_setfield(L, -2, "times");

      // Values table.
      lua_newtable(L);
      for (size_t i = 1; i <= array->size; ++i)
      {
        lua_push_probe_data(L, array->data[i-1]);
        lua_rawseti(L, -2, i);
      }
      lua_setfield(L, -2, "values");

      lua_setfield(L, -2, quantity);
    }
  }
  else
  {
    if (!lua_isstring(L, 2))
      luaL_error(L, "Argument must be the name of probe-acquired data.");
    const char* data_name = lua_tostring(L, 2);
    probe_data_array_t* array = model_probe_data(m, data_name);
    if (array == NULL)
    {
      lua_pushnil(L);
      return 0;
    }

    // We return a table with two fields: times and values.
    lua_newtable(L);

    // Times table.
    lua_newtable(L);
    for (size_t i = 1; i <= array->size; ++i)
    {
      lua_push_real(L, array->data[i-1]->time);
      lua_rawseti(L, -2, i);
    }
    lua_setfield(L, -2, "times");

    // Values table.
    lua_newtable(L);
    for (size_t i = 1; i <= array->size; ++i)
    {
      lua_push_probe_data(L, array->data[i-1]);
      lua_rawseti(L, -2, i);
    }
    lua_setfield(L, -2, "values");
  }
  return 1;
}

static int m_tostring(lua_State* L)
{
  model_t* m = lua_to_model(L, 1);
  if (m == NULL)
    return luaL_error(L, "Method must be invoked with a model.");
  lua_pushfstring(L, "model '%s'", model_name(m));
  return 1;
}

static lua_class_method model_methods[] = {
  {"init", m_init, "model:init(t) -> Initializes the model at time t."},
  {"load", m_load, "model:load(step) -> Loads the given step from a save file."},
  {"save", m_save, "model:save() -> Writes a save file for the given simulation step."},
  {"plot", m_plot, "model:plot() -> Generates a plot for the given simulation step."},
  {"advance", m_advance, "model:advance(max_dt) -> Advances one step, returning the actual timestep taken."},
  {"run", m_run, "model:run{t1 = T1, t2 = T2, max_steps = N, prefix = PREFIX, directory = DIR, plot_every = PLOT_EVERY, save_every = SAVE_EVERY, load_step = LOAD_STEP, min_dt = MIN_DT, max_dt = MAX_DT} -> Runs the model from T1 to T2 or for MAX_STEPS."},
  {"handle_signals", m_handle_signals, "model:handle_signals(flag) -> Sets whether the model intercepts SIGINT and SIGTERM during time stepping. On by default."},
  {"add_probe", m_add_probe, "model:add_probe(probe, times) -> Adds a probe that acquires data at the given list of times."},
  {"data", m_data, "model:data(data_name) -> Returns a table of acquired probe data, or a table of all data if no name is given."},
  {"__tostring", m_tostring, NULL},
  {NULL, NULL, NULL}
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

  // Call the thing with t as the argument.
  ASSERT(lua_istable(p->L, -1) || lua_isfunction(p->L, -1));
  lua_push_real(p->L, t);
  lua_call(p->L, 1, 1);

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
  {"new", p_new, "probe.new(quantity_name, callable_func, data_shape) -> New probe that acquires the given quantity using the given callable function."},
  {NULL, NULL, NULL}
};

static int p_get_name(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  lua_pushfstring(L, probe_name(p));
  return 1;
}

static int p_get_data_name(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  lua_pushstring(L, probe_data_name(p));
  return 1;
}

static lua_class_field probe_fields[] = {
  {"name", p_get_name, NULL},
  {"data_name", p_get_data_name, NULL},
  {NULL, NULL, NULL}
};

static int p_acquire(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  if (p == NULL)
    return luaL_error(L, "Method must be invoked with a probe.");
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_isnumber(L, 2))
    return luaL_error(L, "Argument must be a time.");

  // Do the acquisition.
  real_t t = lua_to_real(L, 2);
  probe_data_t* data = probe_acquire(p, t);

  // Return the acquired data as a table.
  lua_newtable(L);
  lua_push_real(L, t);
  lua_setfield(L, -2, "time");
  lua_push_probe_data(L, data);
  probe_data_free(data);
  lua_setfield(L, -2, "data");
  return 1;
}

// This stuff wraps lua callbacks passed to probe_on_acquire.
typedef struct
{
  lua_State* L;
} lua_acq_func_t;

static void lua_on_acquire(void* context, real_t t, probe_data_t* data)
{
  // Access the Lua function.
  lua_acq_func_t* func = context;
  lua_pushlightuserdata(func->L, func);
  lua_gettable(func->L, LUA_REGISTRYINDEX);
  ASSERT(lua_isfunction(func->L, -1));

  // Push the time and probe data onto the stack.
  lua_push_real(func->L, t);
  lua_push_probe_data(func->L, data);

  // Call the function.
  lua_call(func->L, 2, 0);
}

static int p_on_acquire(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  if (p == NULL)
    return luaL_error(L, "Method must be invoked with a probe.");
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_isfunction(L, 2))
    return luaL_error(L, "Argument must be a function that accepts a time and a probe_data.");

  // Stash this function in Lua's registry.
  lua_acq_func_t* f = polymec_malloc(sizeof(lua_acq_func_t));
  f->L = L;
  lua_pushlightuserdata(L, f);   // key   <-- f
  lua_pushvalue(L, 2);           // value <-- callable object
  lua_settable(L, LUA_REGISTRYINDEX);

  // Add the context and acquisition function to the probe.
  probe_on_acquire(p, f, lua_on_acquire, polymec_free);

  return 0;
}

static int p_stream_on_acquire(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  if (p == NULL)
    return luaL_error(L, "Method must be invoked with a probe.");
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_istable(L, 2))
    return luaL_error(L, "Argument must be a table with 'format' and 'destination' (and possibly 'port') keys.");

  const char* destination = NULL;
  lua_getfield(L, 2, "destination");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "destination must be an URL or a file path.");
  destination = lua_tostring(L, -1);

  int port = -1;
  lua_getfield(L, 2, "port");
  if (!lua_isnil(L, -1))
  {
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "If given, port must be a positive integer.");
    else
    {
      port = (int)lua_tointeger(L, -1);
      if (port <= 0)
        return luaL_error(L, "port must be positive.");
    }
  }

  // If our destination includes a port, use it.
  char host[strlen(destination)+1];
  if (string_num_tokens(destination, ":") == 2)
  {
    int pos = 0;
    size_t len;
    char* host_ptr;
    bool found = string_next_token(destination, ":", &pos, &host_ptr, &len);
    ASSERT(found);
    strncpy(host, host_ptr, len);
    host[len] = '\0';
    char* port_str;
    found = string_next_token(destination, ":", &pos, &port_str, &len);
    ASSERT(found);
    if (!string_is_number(port_str))
      return luaL_error(L, "Invalid port: %s", port_str);
    else
    {
      port = atoi(port_str);
      if (port <= 0)
        return luaL_error(L, "port must be positive.");
    }
  }
  else
    strcpy(host, destination);

  const char* format = NULL;
  lua_getfield(L, 2, "format");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "format must be a string.");
  format = lua_tostring(L, -1);

  bool result = probe_stream_on_acquire(p, host, port, format);
  lua_pushboolean(L, result);
  return 1;
}

static int p_tostring(lua_State* L)
{
  probe_t* p = lua_to_probe(L, 1);
  if (p == NULL)
    return luaL_error(L, "Method must be invoked with a probe.");
  lua_pushfstring(L, "probe '%s' (measures %s)", probe_name(p), probe_data_name(p));
  return 1;
}

static lua_class_method probe_methods[] = {
  {"acquire", p_acquire, "probe:acquire(t) -> Acquires and returns probe data at time t."},
  {"on_acquire", p_on_acquire, "probe:on_acquire(f) -> Adds the function f(t, data) to the set of functions called on each acquisition for this probe."},
  {"stream_on_acquire", p_stream_on_acquire, "probe:stream_on_acquire{destination[, port]} or probe:stream_on_acquire(destination) -> Streams data from the probe to the given destination via UDP or a UNIX socket."},
  {"__tostring", p_tostring, NULL},
  {NULL, NULL, NULL}
};

int lua_register_model_modules(lua_State* L)
{
  lua_register_class(L, "model", "A simulation model.", model_funcs, model_fields, model_methods, DTOR(model_free));
  lua_register_class(L, "probe", "A virtual simulation probe that acquires data.", probe_funcs, probe_fields, probe_methods, DTOR(probe_free));

  return 0;
}

void lua_push_model(lua_State* L, model_t* m)
{
  lua_push_object(L, "model", m);
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
  lua_push_object(L, "probe", p);
}

bool lua_is_probe(lua_State* L, int index)
{
  return lua_is_object(L, index, "probe");
}

probe_t* lua_to_probe(lua_State* L, int index)
{
  return (probe_t*)lua_to_object(L, index, "probe");
}

