// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "io/lua_io.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int silo_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with prefix, dir entries.");

  lua_getfield(L, -2, "prefix");
  if (lua_isnil(L, -1))
    return luaL_error(L, "prefix must be specified.");
  const char* prefix = lua_tostring(L, -1);
  lua_pop(L, 1);

  lua_getfield(L, -2, "dir");
  if (lua_isnil(L, -1))
    return luaL_error(L, "dir must be specified.");
  const char* dir = lua_tostring(L, -1);
  lua_pop(L, 1);

  int num_files = 1;
  lua_getfield(L, -2, "num_files");
  if (!lua_isnil(L, -1))
    num_files = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  int step = 0;
  lua_getfield(L, -2, "step");
  if (!lua_isnil(L, -1))
    step = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  real_t time = 0.0;
  lua_getfield(L, -2, "time");
  if (!lua_isnil(L, -1))
    time = (real_t)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  silo_file_t* s = silo_file_new(MPI_COMM_WORLD, prefix, dir, num_files, step, time);
  lua_push_silo_file(L, s);
  return 1;
}

static int silo_open(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with prefix, dir, step entries.");

  lua_getfield(L, -2, "prefix");
  if (lua_isnil(L, -1))
    return luaL_error(L, "prefix must be specified.");
  const char* prefix = lua_tostring(L, -1);
  lua_pop(L, 1);

  lua_getfield(L, -2, "dir");
  if (lua_isnil(L, -1))
    return luaL_error(L, "dir must be specified.");
  const char* dir = lua_tostring(L, -1);
  lua_pop(L, 1);

  int num_files = 1;
  lua_getfield(L, -2, "num_files");
  if (!lua_isnil(L, -1))
    num_files = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  int step = 0;
  lua_getfield(L, -2, "step");
  if (lua_isnil(L, -1))
    return luaL_error(L, "step must be specified.");
  step = (int)(lua_tonumber(L, -1));
  lua_pop(L, 1);

  real_t time;
  silo_file_t* s = silo_file_open(MPI_COMM_WORLD, prefix, dir, step, &time);
  lua_push_silo_file(L, s);
  return 1;
}

static lua_module_function silo_funcs[] = {
  {"new", silo_new, "silo_file.new{prefix = PREFIX, dir = DIR, num_files = N, step = STEP, time = TIME} -> Opens a SILO file for writing."},
  {"open", silo_open, "silo_file.open{prefix = PREFIX, dir = DIR, num_files = N, step = STEP} -> Opens a SILO file for reading."},
  {NULL, NULL, NULL}
};

static int silo_close(lua_State* L)
{
  silo_file_t* s = lua_to_silo_file(L, 1);
  if (s == NULL)
    luaL_error(L, "Method must be invoked with a silo_file.");
  silo_file_close(s);
  return 0;
}

static lua_class_method silo_methods[] = {
  {"close", silo_close, "file:close() -> Closes the SILO file."},
  {NULL, NULL, NULL}
};

//------------------------------------------------------------------------
//                                API
//------------------------------------------------------------------------

int lua_register_io_modules(lua_State* L)
{
  lua_register_class(L, "silo_file", "A SILO file interface.", silo_funcs, NULL, silo_methods, DTOR(silo_file_close));
  return 0;
}

void lua_push_silo_file(lua_State* L, silo_file_t* s)
{
  lua_push_object(L, "silo_file", s);
}

bool lua_is_silo_file(lua_State* L, int index)
{
  return lua_is_object(L, index, "silo_file");
}

silo_file_t* lua_to_silo_file(lua_State* L, int index)
{
  return (silo_file_t*)lua_to_object(L, index, "silo_file");
}

