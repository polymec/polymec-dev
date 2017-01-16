// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/lua_core.h"
#include "core/partition_mesh.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int p_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 3) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2) || !lua_isnumber(L, 3))
  {
    return luaL_error(L, "Arguments must be x, y, z coordinates.");
  }

  real_t x = (real_t)lua_tonumber(L, 1);
  real_t y = (real_t)lua_tonumber(L, 2);
  real_t z = (real_t)lua_tonumber(L, 3);
  point_t* point = point_new(x, y, z);
  lua_push_point(L, point);
  return 1;
}

static lua_type_func point_funcs[] = {
  {"new", p_new},
  {NULL, NULL}
};

static int p_x(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushnumber(L, (double)p->x);
  return 1;
}

static int p_set_x(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->x = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int p_y(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushnumber(L, (double)p->y);
  return 1;
}

static int p_set_y(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->y = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int p_z(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushnumber(L, (double)p->z);
  return 1;
}

static int p_set_z(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->z = (real_t)lua_tonumber(L, 2);
  return 0;
}

static lua_type_attr point_attr[] = {
  {"x", p_x, p_set_x},
  {"y", p_y, p_set_y},
  {"z", p_z, p_set_z},
  {NULL, NULL, NULL}
};

static int p_distance(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  point_t* q = lua_to_point(L, 2);
  if (q == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushnumber(L, (double)point_distance(p, q));
  return 1;
}

static lua_type_method point_methods[] = {
  {"distance", p_distance},
  {NULL, NULL}
};

static int v_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 3) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2) || !lua_isnumber(L, 3))
  {
    return luaL_error(L, "Arguments must be x, y, z components.");
  }

  real_t x = (real_t)lua_tonumber(L, 1);
  real_t y = (real_t)lua_tonumber(L, 2);
  real_t z = (real_t)lua_tonumber(L, 3);
  vector_t* vec = vector_new(x, y, z);
  lua_push_vector(L, vec);
  return 1;
}

static lua_type_func vector_funcs[] = {
  {"new", v_new},
  {NULL, NULL}
};

static int v_x(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)v->x);
  return 1;
}

static int v_set_x(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Vector components must be numbers.");
  v->x = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int v_y(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)v->y);
  return 1;
}

static int v_set_y(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Vector components must be numbers.");
  v->y = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int v_z(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)v->z);
  return 1;
}

static int v_set_z(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Vector components must be numbers.");
  v->z = (real_t)lua_tonumber(L, 2);
  return 0;
}

static lua_type_attr vector_attr[] = {
  {"x", v_x, v_set_x},
  {"y", v_y, v_set_y},
  {"z", v_z, v_set_z},
  {NULL, NULL, NULL}
};

static int v_dot(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (w == NULL)
    return luaL_error(L, "Argument must be a vector.");
  lua_pushnumber(L, (double)vector_dot(v, w));
  return 1;
}

static lua_type_method vector_methods[] = {
  {"dot", v_dot},
  {NULL, NULL}
};

static int bb_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if (num_args != 1)
  {
    if (!lua_istable(L, 1))
      return luaL_error(L, "Argument must be a table containing x1, x2, y1, y2, z1, z2 values.");
  }

  // Look for x1, x2, y1, y2, z1, z2 in the table.
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  const char* entries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  for (int i = 0; i < 6; ++i)
  {
    lua_pushstring(L, entries[i]);
    lua_gettable(L, 1); // Reads name from top, replaces with bounds[name].
    if (!lua_isnumber(L, -1))
    {
      return luaL_error(L, "Invalid entry for '%s'.\n"
                        "x1, x2, y1, y2, z1, z2, must all be numbers.", entries[i]);
    }
    switch(i)
    {
      case 0: bbox->x1 = (real_t)lua_tonumber(L, -1);
              break;
      case 1: bbox->x2 = (real_t)lua_tonumber(L, -1);
              break;
      case 2: bbox->y1 = (real_t)lua_tonumber(L, -1);
              break;
      case 3: bbox->y2 = (real_t)lua_tonumber(L, -1);
              break;
      case 4: bbox->z1 = (real_t)lua_tonumber(L, -1);
              break;
      case 5: bbox->z2 = (real_t)lua_tonumber(L, -1);
              break;
      default: break;
    }
    lua_pop(L, 1); 
  }

  // Push the bounding box onto the stack.
  lua_push_bbox(L, bbox);
  return 1;
}

static lua_type_func bbox_funcs[] = {
  {"new", bb_new},
  {NULL, NULL}
};

static int bb_x1(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)b->x1);
  return 1;
}

static int bb_set_x1(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding coordinates must be numbers.");
  real_t x1 = (real_t)lua_tonumber(L, 2);
  if (x1 > b->x2) 
    return luaL_error(L, "Bounding intervals must be positive.");
  b->x1 = x1;
  return 0;
}

static int bb_x2(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)b->x2);
  return 1;
}

static int bb_set_x2(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding coordinates must be numbers.");
  real_t x2 = (real_t)lua_tonumber(L, 2);
  if (x2 < b->x1) 
    return luaL_error(L, "Bounding intervals must be positive.");
  b->x2 = x2;
  return 0;
}

static int bb_y1(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)b->y1);
  return 1;
}

static int bb_set_y1(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding coordinates must be numbers.");
  real_t y1 = (real_t)lua_tonumber(L, 2);
  if (y1 > b->y2) 
    return luaL_error(L, "Bounding intervals must be positive.");
  b->y1 = y1;
  return 0;
}

static int bb_y2(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)b->y2);
  return 1;
}

static int bb_set_y2(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding coordinates must be numbers.");
  real_t y2 = (real_t)lua_tonumber(L, 2);
  if (y2 < b->y1) 
    return luaL_error(L, "Bounding intervals must be positive.");
  b->y2 = y2;
  return 0;
}

static int bb_z1(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)b->z1);
  return 1;
}

static int bb_set_z1(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding coordinates must be numbers.");
  real_t z1 = (real_t)lua_tonumber(L, 2);
  if (z1 > b->z2) 
    return luaL_error(L, "Bounding intervals must be positive.");
  b->z1 = z1;
  return 0;
}

static int bb_z2(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)b->z2);
  return 1;
}

static int bb_set_z2(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding coordinates must be numbers.");
  real_t z2 = (real_t)lua_tonumber(L, 2);
  if (z2 < b->z1) 
    return luaL_error(L, "Bounding intervals must be positive.");
  b->z2 = z2;
  return 0;
}

static lua_type_attr bbox_attr[] = {
  {"x1", bb_x1, bb_set_x1},
  {"x2", bb_x2, bb_set_x2},
  {"x1", bb_y1, bb_set_y1},
  {"x2", bb_y2, bb_set_y2},
  {"x1", bb_z1, bb_set_z1},
  {"x2", bb_z2, bb_set_z2},
  {NULL, NULL, NULL}
};

static int bb_contains(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  point_t* p = lua_to_point(L, 2);
  if (p == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushboolean(L, bbox_contains(b, p));
  return 1;
}

static lua_type_method bbox_methods[] = {
  {"contains", bb_contains},
  {NULL, NULL}
};

static int sp_constant(lua_State* L)
{
  // Check the argument.
  int num_args = lua_gettop(L);
  real_t val[num_args];
  for (int i = 0; i < num_args; ++i)
  {
    if (!lua_isnumber(L, i))
      return luaL_error(L, "Argument %d must be a number.", i);
    val[i] = (real_t)lua_tonumber(L, i);
  }
  sp_func_t* f = constant_sp_func_new(val, num_args);
  lua_push_sp_func(L, f);
  return 1;
}

static lua_type_func sp_funcs[] = {
  {"constant", sp_constant},
  {NULL, NULL}
};

static int sp_num_comp(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  lua_pushnumber(L, (double)sp_func_num_comp(f));
  return 1;
}

static lua_type_attr sp_attr[] = {
  {"num_comp", sp_num_comp, NULL},
  {NULL, NULL, NULL}
};

static int sp_call(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  int nc = sp_func_num_comp(f);
  real_t val[nc];
  sp_func_eval(f, x, val);
  for (int i = 0; i < nc; ++i)
    lua_pushnumber(L, val[i]);
  return nc;
}

static lua_type_method sp_methods[] = {
  {"__call", sp_call},
  {NULL, NULL}
};

static int st_constant(lua_State* L)
{
  // Check the argument.
  int num_args = lua_gettop(L);
  real_t val[num_args];
  for (int i = 0; i < num_args; ++i)
  {
    if (!lua_isnumber(L, i))
      return luaL_error(L, "Argument %d must be a number.", i);
    val[i] = (real_t)lua_tonumber(L, i);
  }
  st_func_t* f = constant_st_func_new(val, num_args);
  lua_push_st_func(L, f);
  return 1;
}

static lua_type_func st_funcs[] = {
  {"constant", st_constant},
  {NULL, NULL}
};

static int st_num_comp(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  lua_pushnumber(L, (double)st_func_num_comp(f));
  return 1;
}

static lua_type_attr st_attr[] = {
  {"num_comp", st_num_comp, NULL},
  {NULL, NULL, NULL}
};

static int st_call(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "First argument must be a point.");
  if (!lua_is_point(L, 3))
    return luaL_error(L, "Second argument must be a time.");
  point_t* x = lua_to_point(L, 2);
  real_t t = (real_t)lua_tonumber(L, 3);
  int nc = st_func_num_comp(f);
  real_t val[nc];
  st_func_eval(f, x, t, val);
  for (int i = 0; i < nc; ++i)
    lua_pushnumber(L, val[i]);
  return nc;
}

static lua_type_method st_methods[] = {
  {"__call", st_call},
  {NULL, NULL}
};

static int mesh_repartition(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_mesh(L, 1))
    return luaL_error(L, "Argument must be mesh.");
  mesh_t* mesh = lua_to_mesh(L, 1);
  real_t imbalance_tol = 0.05;

  // Bug out if there's only one process.
  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return 0;

  // Make sure there are enough cells for our processes.
  index_t local_num_cells = mesh->num_cells, global_num_cells;
  MPI_Allreduce(&local_num_cells, &global_num_cells, 1, MPI_INDEX_T, MPI_SUM, mesh->comm);
  if (global_num_cells < nprocs)
    return luaL_error(L, "Insufficient number of cells (%zd) for number of processes (%d).", global_num_cells, nprocs);

  // Perform the repartitioning and toss the migrator.
  migrator_t* m = repartition_mesh(&mesh, NULL, imbalance_tol);
  m = NULL;

  return 0;
}

static lua_type_func mesh_funcs[] = {
  {"repartition", mesh_repartition},
  {NULL, NULL}
};

static int mesh_num_cells(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushnumber(L, (double)m->num_cells);
  return 1;
}

static lua_type_attr mesh_attr[] = {
  {"num_cells", mesh_num_cells, NULL},
  {NULL, NULL, NULL}
};

static lua_type_method mesh_methods[] = {
  {NULL, NULL}
};

static int pc_repartition(lua_State* L)
{
  luaL_error(L, "can't repartion point clouds just yet!");
  return 0;
}

static lua_type_func pc_funcs[] = {
  {"repartition", pc_repartition},
  {NULL, NULL}
};

static int pc_num_points(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushnumber(L, (double)pc->num_points);
  return 1;
}

static lua_type_attr pc_attr[] = {
  {"num_points", pc_num_points, NULL},
  {NULL, NULL, NULL}
};

static lua_type_method pc_methods[] = {
  {NULL, NULL}
};

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
    num_files = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);

  int step = 0;
  lua_getfield(L, -2, "step");
  if (!lua_isnil(L, -1))
    step = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);

  real_t time = 0.0;
  lua_getfield(L, -2, "time");
  if (!lua_isnil(L, -1))
    time = (real_t)lua_tonumber(L, -1);
  lua_pop(L, 1);

  silo_file_t* s = silo_file_new(MPI_COMM_WORLD, prefix, dir, num_files, 0, step, time);
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
    num_files = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);

  int step = 0;
  lua_getfield(L, -2, "step");
  if (lua_isnil(L, -1))
    return luaL_error(L, "step must be specified.");
  step = (int)lua_tonumber(L, -1);
  lua_pop(L, 1);

  real_t time;
  silo_file_t* s = silo_file_open(MPI_COMM_WORLD, prefix, dir, 0, step, &time);
  lua_push_silo_file(L, s);
  return 1;
}

static lua_type_func silo_funcs[] = {
  {"new", silo_new},
  {"open", silo_open},
  {NULL, NULL}
};


static lua_type_attr silo_attr[] = {
  {NULL, NULL, NULL}
};

static int silo_close(lua_State* L)
{
  silo_file_t* s = lua_to_silo_file(L, 1);
  silo_file_close(s);
  return 0;
}

static lua_type_method silo_methods[] = {
  {"close", silo_close},
  {NULL, NULL}
};

void lua_register_core_modules(lua_State* L)
{
  // Core types.
  lua_register_type(L, "point", point_funcs, point_attr, point_methods);
  lua_register_type(L, "vector", vector_funcs, vector_attr, vector_methods);
  lua_register_type(L, "bbox", bbox_funcs, bbox_attr, bbox_methods);
  lua_register_type(L, "sp_func", sp_funcs, sp_attr, sp_methods);
  lua_register_type(L, "st_func", st_funcs, st_attr, st_methods);
  lua_register_type(L, "mesh", mesh_funcs, mesh_attr, mesh_methods);
  lua_register_type(L, "point_cloud", pc_funcs, pc_attr, pc_methods);
  lua_register_type(L, "silo_file", silo_funcs, silo_attr, silo_methods);
}

void lua_push_point(lua_State* L, point_t* p)
{
  lua_push_object(L, "point", p, NULL);
}

bool lua_is_point(lua_State* L, int index)
{
  return lua_is_object(L, index, "point");
}

point_t* lua_to_point(lua_State* L, int index)
{
  return (point_t*)lua_to_object(L, index, "point");
}

void lua_push_vector(lua_State* L, vector_t* v)
{
  lua_push_object(L, "vector", v, NULL);
}

bool lua_is_vector(lua_State* L, int index)
{
  return lua_is_object(L, index, "vector");
}

vector_t* lua_to_vector(lua_State* L, int index)
{
  return (vector_t*)lua_to_object(L, index, "vector");
}

void lua_push_bbox(lua_State* L, bbox_t* b)
{
  lua_push_object(L, "bbox", b, NULL);
}

bool lua_is_bbox(lua_State* L, int index)
{
  return lua_is_object(L, index, "bbox");
}

bbox_t* lua_to_bbox(lua_State* L, int index)
{
  return (bbox_t*)lua_to_object(L, index, "bbox");
}

void lua_push_sp_func(lua_State* L, sp_func_t* f)
{
  lua_push_object(L, "sp_func", f, NULL);
}

bool lua_is_sp_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sp_func");
}

sp_func_t* lua_to_sp_func(lua_State* L, int index)
{
  return (sp_func_t*)lua_to_object(L, index, "sp_func");
}

void lua_push_st_func(lua_State* L, st_func_t* f)
{
  lua_push_object(L, "st_func", f, NULL);
}

bool lua_is_st_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "st_func");
}

st_func_t* lua_to_st_func(lua_State* L, int index)
{
  return (st_func_t*)lua_to_object(L, index, "st_func");
}

void lua_push_mesh(lua_State* L, mesh_t* m)
{
  lua_push_object(L, "mesh", m, DTOR(mesh_free));
}

bool lua_is_mesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "mesh");
}

mesh_t* lua_to_mesh(lua_State* L, int index)
{
  return (mesh_t*)lua_to_object(L, index, "mesh");
}

void lua_push_point_cloud(lua_State* L, point_cloud_t* c)
{
  lua_push_object(L, "point_cloud", c, DTOR(point_cloud_free));
}

bool lua_is_point_cloud(lua_State* L, int index)
{
  return lua_is_object(L, index, "point_cloud");
}

point_cloud_t* lua_to_point_cloud(lua_State* L, int index)
{
  return (point_cloud_t*)lua_to_object(L, index, "point_cloud");
}

void lua_push_silo_file(lua_State* L, silo_file_t* s)
{
  lua_push_object(L, "silo_file", s, DTOR(silo_file_close));
}

bool lua_is_silo_file(lua_State* L, int index)
{
  return lua_is_object(L, index, "silo_file");
}

silo_file_t* lua_to_silo_file(lua_State* L, int index)
{
  return (silo_file_t*)lua_to_object(L, index, "silo_file");
}
