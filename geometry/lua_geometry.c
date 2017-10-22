// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"
#include "geometry/lua_geometry.h"
#include "geometry/plane_sd_func.h"
#include "geometry/sphere_sd_func.h"
#include "geometry/cylinder_sd_func.h"
#include "geometry/union_sd_func.h"
#include "geometry/intersection_sd_func.h"
#include "geometry/difference_sd_func.h"

#include "geometry/partition_polymesh.h"
#include "geometry/create_uniform_polymesh.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This serves as a context pointer that houses Lua objects.
typedef struct
{
  lua_State* L;
} lua_obj_t;

static void cm_map_point(void* context, point_t* x, point_t* y)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the map_point method.
  lua_getfield(lo->L, -1, "map_point");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_call(lo->L, 2, 1);

  if (!lua_is_point(lo->L, -1))
    luaL_error(lo->L, "map_point method did not return a point.");
  *y = *(lua_to_point(lo->L, -1));
}

static void cm_compute_jacobian(void* context, point_t* x, tensor2_t* J)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the jacobian method.
  lua_getfield(lo->L, -1, "jacobian");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_call(lo->L, 2, 1);

  if (!lua_is_tensor2(lo->L, -1))
    luaL_error(lo->L, "jacobian method did not return a tensor2.");
  *J = *(lua_to_tensor2(lo->L, -1));
}

static coord_mapping_t* cm_compute_inverse(void* context)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the jacobian method.
  lua_getfield(lo->L, -1, "inverse");
  lua_pushvalue(lo->L, -2);
  lua_call(lo->L, 1, 1);

  if (!lua_is_coord_mapping(lo->L, -1))
  {
    luaL_error(lo->L, "inverse method did not return a coord_mapping.");
    return NULL; 
  }
  return lua_to_coord_mapping(lo->L, -1);
}

static void cm_dtor(void* context)
{
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);
  lua_pushnil(lo->L);
  lua_settable(lo->L, LUA_REGISTRYINDEX);
  polymec_free(context);
}

static int cm_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with a name and methods.");

  // Create a context that houses this object.
  lua_newtable(L);
  int obj_index = num_args + 1;

  coord_mapping_vtable vtable = {.dtor = cm_dtor};

  lua_getfield(L, 1, "name");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "Name must be a string.");
  const char* name = lua_tostring(L, -1);

  lua_getfield(L, 1, "map");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "map must be a method.");
  lua_setfield(L, obj_index, "map");
  vtable.map_point = cm_map_point;

  lua_getfield(L, 1, "jacobian");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "jacobian must be a method.");
  lua_setfield(L, obj_index, "jacobian");
  vtable.jacobian = cm_compute_jacobian;

  lua_getfield(L, 1, "inverse");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "jacobian must be a method.");
  lua_setfield(L, obj_index, "jacobian");
  vtable.inverse = cm_compute_inverse;

  // Allocate a context pointer and stuff our object into the registry.
  lua_obj_t* lo = polymec_malloc(sizeof(lua_obj_t));
  lo->L = L;

  // Store the table representing our object in the registry, with 
  // context as a key.
  lua_pushvalue(L, obj_index);
  lua_rawsetp(L, LUA_REGISTRYINDEX, lo);

  // Set up the mapping.
  coord_mapping_t* X = coord_mapping_new(name, lo, vtable);
  lua_push_coord_mapping(L, X);

  return 1;
}

static int cm_compose(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_is_coord_mapping(L, 2) || !lua_is_coord_mapping(L, 3))
    return luaL_error(L, "Arguments must both be coord_mappings.");

  coord_mapping_t* X1 = lua_to_coord_mapping(L, 2);
  coord_mapping_t* X2 = lua_to_coord_mapping(L, 3);
  coord_mapping_t* X1oX2 = composite_coord_mapping_new(X1, X2);
  lua_push_coord_mapping(L, X1oX2);
  return 1;
}

static lua_module_function cm_funcs[] = {
  {"new", cm_new, "coord_mapping.new{name = NAME, map = Xp, jacobian = J, inverse = Xinv} -> new coordinate mapping."},
  {"compose", cm_compose, "coord_mapping.compose(X1, X2) ­> Returns a coordinate mapping for X1 o X2."},
  {NULL, NULL, NULL}
};

static int cm_jacobian(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  if (!lua_is_point(L, 2))
    luaL_error(L, "Argument must be a point.");

  point_t* x = lua_to_point(L, 2);
  tensor2_t* J = tensor2_new(0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0);
  coord_mapping_compute_jacobian(X, x, J);
  lua_push_tensor2(L, J);
  return 1;
}

static int cm_metric(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  if (!lua_is_point(L, 2))
    luaL_error(L, "Argument must be a point.");

  point_t* x = lua_to_point(L, 2);
  tensor2_t* G = tensor2_new(0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0);
  coord_mapping_compute_metric(X, x, G);
  lua_push_tensor2(L, G);
  return 1;
}

static int cm_inverse(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  coord_mapping_t* Xinv = coord_mapping_inverse(X);
  if (Xinv != NULL)
    lua_push_coord_mapping(L, Xinv);
  else
    lua_pushnil(L);
  return 1;
}

static int cm_call(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  int num_args = lua_gettop(L);
  if (!((num_args == 2) && lua_is_point(L, 2)) ||
       ((num_args == 3) && lua_is_point(L, 2) && lua_is_vector(L, 3)))
    luaL_error(L, "Argument must be a point or a point and a vector.");

  if (lua_is_point(L, 2))
  {
    point_t* x = lua_to_point(L, 2);
    point_t* y = point_new(0.0, 0.0, 0.0);
    coord_mapping_map_point(X, x, y);
    lua_push_point(L, y);
  }
  else // vector
  {
    point_t* x = lua_to_point(L, 2);
    vector_t* v = lua_to_vector(L, 3);
    vector_t* v1 = vector_new(0.0, 0.0, 0.0);
    coord_mapping_map_vector(X, x, v, v1);
    lua_push_vector(L, v1);
  }

  return 1;
}

static int cm_tostring(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  lua_pushfstring(L, "coord_mapping '%s'", coord_mapping_name(X));
  return 1;
}

static lua_class_method cm_methods[] = {
  {"jacobian", cm_jacobian, "X:jacobian(x) -> Returns the jacobian of X at the point x."},
  {"inverse", cm_inverse, "X:inverse() -> Returns the inverse mapping of X."},
  {"metric", cm_metric, "X:metric(x) -> Returns the 3x3 metric tensor for X at the point x."},
  {"__call", cm_call, NULL},
  {"__tostring", cm_tostring, NULL},
  {NULL, NULL, NULL}
};

static real_t sd_value(void* context, point_t* x)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the value method.
  lua_getfield(lo->L, -1, "value");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_call(lo->L, 2, 1);

  if (!lua_isnumber(lo->L, -1))
    return luaL_error(lo->L, "value method did not return a number.");
  return lua_to_real(lo->L, -1);
}

static void sd_eval_grad(void* context, point_t* x, vector_t* grad)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the grad method.
  lua_getfield(lo->L, -1, "grad");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_call(lo->L, 2, 1);

  if (!lua_is_vector(lo->L, -1))
    luaL_error(lo->L, "grad method did not return a vector.");
  *grad = *lua_to_vector(lo->L, -1);
}

static void sd_dtor(void* context)
{
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);
  lua_pushnil(lo->L);
  lua_settable(lo->L, LUA_REGISTRYINDEX);
  polymec_free(context);
}

static int sd_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with a name and methods.");

  // Create a context that houses this object.
  lua_newtable(L);
  int obj_index = num_args + 1;

  sd_func_vtable vtable = {.dtor = sd_dtor};

  lua_getfield(L, 1, "name");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "Name must be a string.");
  const char* name = lua_tostring(L, -1);

  lua_getfield(L, 1, "value");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "value must be a method.");
  lua_setfield(L, obj_index, "value");
  vtable.value = sd_value;

  lua_getfield(L, 1, "grad");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "grad must be a method.");
  lua_setfield(L, obj_index, "grad");
  vtable.eval_grad = sd_eval_grad;

  // Allocate a context pointer and stuff our object into the registry.
  lua_obj_t* lo = polymec_malloc(sizeof(lua_obj_t));
  lo->L = L;

  // Store the table representing our object in the registry, with 
  // context as a key.
  lua_pushvalue(L, obj_index);
  lua_rawsetp(L, LUA_REGISTRYINDEX, lo);

  // Set up the function.
  sd_func_t* f = sd_func_new(name, lo, vtable);
  lua_push_sd_func(L, f);

  return 1;
}

static int sd_from_sp_funcs(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    return luaL_error(L, "Arguments must be a name, a signed distance function and its gradient.");
  if (!lua_isstring(L, 1))
    return luaL_error(L, "Argument 1 must be a name.");
  if (!lua_is_sp_func(L, 2))
    return luaL_error(L, "Argument 2 must be a signed distance function.");
  if (!lua_is_sp_func(L, 3))
    return luaL_error(L, "Argument 3 must be the gradient of argument 2.");
  const char* name = lua_tostring(L, 1);
  sp_func_t* D = lua_to_sp_func(L, 2);
  sp_func_t* G = lua_to_sp_func(L, 3);
  lua_push_sd_func(L, sd_func_from_sp_funcs(name, D, G));
  return 1;
}

static int sd_plane(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with normal and point entries.");

  lua_getfield(L, 1, "normal");
  if (!lua_is_vector(L, -1))
    return luaL_error(L, "normal must be a vector.");
  vector_t* n = lua_to_vector(L, -1);

  lua_getfield(L, 1, "point");
  if (!lua_is_point(L, -1))
    return luaL_error(L, "point must be a point.");
  point_t* x = lua_to_point(L, -1);

  sd_func_t* f = plane_sd_func_new(n, x);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_sphere(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with point, radius, and orientation entries.");

  lua_getfield(L, 1, "point");
  if (!lua_is_point(L, -1))
    return luaL_error(L, "point must be a point.");
  point_t* x = lua_to_point(L, -1);

  lua_getfield(L, 1, "radius");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "radius must be a positive number.");
  real_t r = lua_to_real(L, -1);
  if (r <= 0.0)
    return luaL_error(L, "radius must be positive.");

  normal_orient_t orient;
  lua_getfield(L, 1, "orientation");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "orientation must be 'inward' or 'outward'.");
  const char* orientation = lua_tostring(L, -1);
  if (strcasecmp(orientation, "inward") == 0)
    orient = INWARD_NORMAL;
  else if (strcasecmp(orientation, "outward") == 0)
    orient = OUTWARD_NORMAL;
  else 
    return luaL_error(L, "invalid orientation: %s (must be 'inward' or 'outward').", orientation);

  sd_func_t* f = sphere_sd_func_new(x, r, orient);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_cylinder(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with point, radius, and orientation entries.");

  lua_getfield(L, 1, "point");
  if (!lua_is_point(L, -1))
    return luaL_error(L, "point must be a point.");
  point_t* x = lua_to_point(L, -1);

  lua_getfield(L, 1, "radius");
  if (!lua_isnumber(L, -1))
    return luaL_error(L, "radius must be a positive number.");
  real_t r = lua_to_real(L, -1);
  if (r <= 0.0)
    return luaL_error(L, "radius must be positive.");

  normal_orient_t orient;
  lua_getfield(L, 1, "orientation");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "orientation must be 'inward' or 'outward'.");
  const char* orientation = lua_tostring(L, -1);
  if (strcasecmp(orientation, "inward") == 0)
    orient = INWARD_NORMAL;
  else if (strcasecmp(orientation, "outward") == 0)
    orient = OUTWARD_NORMAL;
  else 
    return luaL_error(L, "invalid orientation: %s (must be 'inward' or 'outward').", orientation);

  sd_func_t* f = cylinder_sd_func_new(x, r, orient);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_union(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table of sd_funcs.");

  int num_surfaces = (int)lua_rawlen(L, 1);
  sd_func_t* surfaces[num_surfaces];

  for (int i = 0; i < num_surfaces; ++i)
  {
    lua_rawgeti(L, 1, i+1);
    if (!lua_is_sd_func(L, -1))
      return luaL_error(L, "Item %d in table is not an sd_func.", i+1);
    surfaces[i] = lua_to_sd_func(L, -1);
  }

  sd_func_t* f = union_sd_func_new(surfaces, num_surfaces);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_intersection(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table of sd_funcs.");

  int num_surfaces = (int)lua_rawlen(L, 1);
  sd_func_t* surfaces[num_surfaces];

  for (int i = 0; i < num_surfaces; ++i)
  {
    lua_rawgeti(L, 1, i+1);
    if (!lua_is_sd_func(L, -1))
      return luaL_error(L, "Item %d in table is not an sd_func.", i+1);
    surfaces[i] = lua_to_sd_func(L, -1);
  }

  sd_func_t* f = intersection_sd_func_new(surfaces, num_surfaces);
  lua_push_sd_func(L, f);
  return 1;
}

static int sd_difference(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2)
    return luaL_error(L, "Arguments must be a pair of sd_funcs.");
  if (!lua_is_sd_func(L, 1))
    return luaL_error(L, "Argument 1 must be an sd_func.");
  if (!lua_is_sd_func(L, 2))
    return luaL_error(L, "Argument 2 must be an sd_func.");

  sd_func_t* f = difference_sd_func_new(lua_to_sd_func(L, 1), 
                                        lua_to_sd_func(L, 2));
  lua_push_sd_func(L, f);
  return 1;
}

static lua_module_function sd_funcs[] = {
  {"new", sd_new, "sd_func.new{name = NAME, value = F, grad = gradF} -> new signed distance function."},
  {"from_sp_funcs", sd_from_sp_funcs, "sd_func.from_sp_funcs(name, D, G) -> new signed distance function using D for distance and G for gradients."},
  {"plane", sd_plane, "sd_func.plane{normal = n, point = x} -> new plane signed distance function."},
  {"sphere", sd_sphere, "sd_func.sphere{point = x, radius = r, orientation = 'inward'|'outward'} -> new sphere signed distance function."},
  {"cylinder", sd_cylinder, "sd_func.cylinder{point = x, radius = r, orientation = 'inward'|'outward'} -> new (infinite) cylinder signed distance function."},
  {"union", sd_union, "sd_func.union(funcs) -> new signed distance function (union of funcs)."},
  {"intersection", sd_intersection, "sd_func.intersection(funcs) -> new signed distance function (intersection of funcs)."},
  {"difference", sd_difference, "sd_func.difference(f, g) -> new signed distance function (f - g)."},
  {NULL, NULL, NULL}
};

static int sd_rename(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sd_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static int sd_grad(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* x = lua_to_point(L, 2);
  vector_t* grad = vector_new(0.0, 0.0, 0.0);
  sd_func_eval_grad(f, x, grad);
  lua_push_vector(L, grad);
  return 1;
}

static int sd_call(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  lua_pushnumber(L, sd_func_value(f, x));
  return 1;
}

static int sd_project(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  point_t* proj_x = point_new(0.0, 0.0, 0.0);
  sd_func_project(f, x, proj_x);
  lua_push_point(L, proj_x);
  return 1;
}

static int sd_tostring(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  lua_pushfstring(L, "sd_func '%s'", sd_func_name(f));
  return 1;
}

static lua_class_method sd_methods[] = {
  {"rename", sd_rename, "f:rename(name) -> Renames f to the given name."},
  {"grad", sd_grad, "f:grad(x) -> Returns the gradient of f at x."},
  {"project", sd_project, "f:project(x) -> Returns the projection of f to x."},
  {"__call", sd_call, NULL},
  {"__tostring", sd_tostring, NULL},
  {NULL, NULL, NULL}
};

static real_t sdt_value(void* context, point_t* x, real_t t)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the value method.
  lua_getfield(lo->L, -1, "value");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_push_real(lo->L, t);
  lua_call(lo->L, 3, 1);

  if (!lua_isnumber(lo->L, -1))
    return luaL_error(lo->L, "value method did not return a number.");
  return lua_to_real(lo->L, -1);
}

static void sdt_eval_grad(void* context, point_t* x, real_t t, vector_t* grad)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the grad method.
  lua_getfield(lo->L, -1, "grad");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_push_real(lo->L, t);
  lua_call(lo->L, 3, 1);

  if (!lua_is_vector(lo->L, -1))
    luaL_error(lo->L, "grad method did not return a vector.");
  *grad = *lua_to_vector(lo->L, -1);
}

static void sdt_dtor(void* context)
{
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);
  lua_pushnil(lo->L);
  lua_settable(lo->L, LUA_REGISTRYINDEX);
  polymec_free(context);
}

static int sdt_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table with a name and methods.");

  // Create a context that houses this object.
  lua_newtable(L);
  int obj_index = num_args + 1;

  sdt_func_vtable vtable = {.dtor = sdt_dtor};

  lua_getfield(L, 1, "name");
  if (!lua_isstring(L, -1))
    return luaL_error(L, "Name must be a string.");
  const char* name = lua_tostring(L, -1);

  lua_getfield(L, 1, "value");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "value must be a method.");
  lua_setfield(L, obj_index, "value");
  vtable.value = sdt_value;

  lua_getfield(L, 1, "grad");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "grad must be a method.");
  lua_setfield(L, obj_index, "grad");
  vtable.eval_grad = sdt_eval_grad;

  // Allocate a context pointer and stuff our object into the registry.
  lua_obj_t* lo = polymec_malloc(sizeof(lua_obj_t));
  lo->L = L;

  // Store the table representing our object in the registry, with 
  // context as a key.
  lua_pushvalue(L, obj_index);
  lua_rawsetp(L, LUA_REGISTRYINDEX, lo);

  // Set up the function.
  sdt_func_t* f = sdt_func_new(name, lo, vtable);
  lua_push_sdt_func(L, f);

  return 1;
}

static int sdt_from_st_funcs(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    return luaL_error(L, "Arguments must be a name, a signed distance function and its gradient.");
  if (!lua_isstring(L, 1))
    return luaL_error(L, "Argument 1 must be a name.");
  if (!lua_is_st_func(L, 2))
    return luaL_error(L, "Argument 2 must be a time-dependent signed distance function.");
  if (!lua_is_st_func(L, 3))
    return luaL_error(L, "Argument 3 must be the gradient of argument 2.");
  const char* name = lua_tostring(L, 1);
  st_func_t* D = lua_to_st_func(L, 2);
  st_func_t* G = lua_to_st_func(L, 3);
  lua_push_sdt_func(L, sdt_func_from_st_funcs(name, D, G));
  return 1;
}

static lua_module_function sdt_funcs[] = {
  {"new", sdt_new, "sdt_func.new{name = NAME, value = F, grad = gradF} -> new signed distance function."},
  {"from_st_funcs", sdt_from_st_funcs, "sdt_func.from_st_funcs(name, D, G) -> new signed distance function using D for distance and G for gradients."},
  {NULL, NULL, NULL}
};

static int sdt_rename(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sdt_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static int sdt_grad(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* x = lua_to_point(L, 2);
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Argument 2 must be a time.");
  real_t t = lua_to_real(L, 3);
  vector_t* grad = vector_new(0.0, 0.0, 0.0);
  sdt_func_eval_grad(f, x, t, grad);
  lua_push_vector(L, grad);
  return 1;
}

static int sdt_call(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* x = lua_to_point(L, 2);
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Argument 2 must be a time.");
  real_t t = lua_to_real(L, 3);
  lua_pushnumber(L, sdt_func_value(f, x, t));
  return 1;
}

static int sdt_project(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument 1 must be a point.");
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Argument 2 must be a time.");
  point_t* x = lua_to_point(L, 2);
  real_t t = lua_to_real(L, 3);
  point_t* proj_x = point_new(0.0, 0.0, 0.0);
  sdt_func_project(f, x, t, proj_x);
  lua_push_point(L, proj_x);
  return 1;
}

static int sdt_tostring(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  lua_pushfstring(L, "sdt_func '%s'", sdt_func_name(f));
  return 1;
}

static lua_class_method sdt_methods[] = {
  {"rename", sdt_rename, "f:rename(name) -> Renames f to the given name."},
  {"grad", sdt_grad, "f:grad(x) -> Returns the gradient of f at x."},
  {"project", sdt_project, "f:project(x) -> Returns the projection of f to x."},
  {"__call", sdt_call, NULL},
  {"__tostring", sdt_tostring, NULL},
  {NULL, NULL, NULL}
};

static int polymesh_repartition(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_polymesh(L, 1))
    return luaL_error(L, "Argument must be a polymesh.");
  polymesh_t* mesh = lua_to_polymesh(L, 1);
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
  migrator_t* m = repartition_polymesh(&mesh, NULL, imbalance_tol);
  m = NULL;

  return 0;
}

static lua_module_function polymesh_funcs[] = {
  {"repartition", polymesh_repartition, "mesh.repartition(m) -> Repartitions the polymesh m."},
  {NULL, NULL, NULL}
};

static int polymesh_num_cells(lua_State* L)
{
  polymesh_t* m = lua_to_polymesh(L, 1);
  lua_pushinteger(L, m->num_cells);
  return 1;
}

static int polymesh_num_ghost_cells(lua_State* L)
{
  polymesh_t* m = lua_to_polymesh(L, 1);
  lua_pushinteger(L, m->num_ghost_cells);
  return 1;
}

static int polymesh_num_faces(lua_State* L)
{
  polymesh_t* m = lua_to_polymesh(L, 1);
  lua_pushinteger(L, m->num_faces);
  return 1;
}

static int polymesh_num_edges(lua_State* L)
{
  polymesh_t* m = lua_to_polymesh(L, 1);
  lua_pushinteger(L, m->num_edges);
  return 1;
}

static int polymesh_num_nodes(lua_State* L)
{
  polymesh_t* m = lua_to_polymesh(L, 1);
  lua_pushinteger(L, m->num_nodes);
  return 1;
}

static lua_record_field polymesh_fields[] = {
  {"num_cells", polymesh_num_cells, NULL},
  {"num_ghost_cells", polymesh_num_ghost_cells, NULL},
  {"num_faces", polymesh_num_faces, NULL},
  {"num_edges", polymesh_num_edges, NULL},
  {"num_nodes", polymesh_num_nodes, NULL},
  {NULL, NULL, NULL}
};

static int polymesh_tostring(lua_State* L)
{
  polymesh_t* m = lua_to_polymesh(L, 1);
  lua_pushfstring(L, "polymesh (%d cells, %d faces, %d nodes)", 
                  m->num_cells, m->num_faces, m->num_nodes);
  return 1;
}

static lua_record_metamethod polymesh_mm[] = {
  {"__tostring", polymesh_tostring},
  {NULL, NULL}
};

static int polymeshes_uniform(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, nx, ny, nz, bbox fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm object.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "nx");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "nx must be a positive integer.");
  int nx = (int)lua_tointeger(L, -1);
  if (nx < 1)
    luaL_error(L, "nx must be positive.");

  lua_getfield(L, 1, "ny");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "ny must be a positive integer.");
  int ny = (int)lua_tointeger(L, -1);
  if (ny < 1)
    luaL_error(L, "ny must be positive.");

  lua_getfield(L, 1, "nz");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "nz must be a positive integer.");
  int nz = (int)lua_tointeger(L, -1);
  if (nz < 1)
    luaL_error(L, "nz must be positive.");

  lua_getfield(L, 1, "bbox");
  if (!lua_is_bbox(L, -1))
    luaL_error(L, "bbox must be a bounding box (bbox).");
  bbox_t* bbox = lua_to_bbox(L, -1);

  int rank = -1;
  lua_getfield(L, 1, "rank");
  if (lua_isinteger(L, -1))
    rank = (int)lua_tointeger(L, -1);

  polymesh_t* mesh = NULL;
  if (rank == -1) 
    mesh = create_uniform_polymesh(comm, nx, ny, nz, bbox);
  else
    mesh = create_uniform_polymesh_on_rank(comm, rank, nx, ny, nz, bbox);

  lua_push_polymesh(L, mesh);
  return 1;
}

static real_array_t* get_coordinates(lua_State* L, int index)
{
  real_array_t* array = NULL;
  if (lua_is_array(L, index, LUA_ARRAY_REAL))
    array = lua_to_array(L, index, LUA_ARRAY_REAL);
  else
  {
    array = real_array_new();
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, index, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_is_real(L, -1))
          luaL_error(L, "Item %d in table is not a coordinate.");
        real_array_append(array, lua_to_real(L, -1));
        lua_pop(L, 1);
      }
      ++i;
    }
  }
  return array;
}

static int polymeshes_rectilinear(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, xs, ys, zs fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm object.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "xs");
  if (!lua_is_array(L, -1, LUA_ARRAY_REAL) && !lua_istable(L, -1))
    luaL_error(L, "xs must be a table or array of x coordinates.");
  real_array_t* xs = get_coordinates(L, -1);

  lua_getfield(L, 1, "ys");
  if (!lua_is_array(L, -1, LUA_ARRAY_REAL) && !lua_istable(L, -1))
    luaL_error(L, "ys must be a table or array of y coordinates.");
  real_array_t* ys = get_coordinates(L, -1);

  lua_getfield(L, 1, "zs");
  if (!lua_is_array(L, -1, LUA_ARRAY_REAL) && !lua_istable(L, -1))
    luaL_error(L, "zs must be a table or array of z coordinates.");
  real_array_t* zs = get_coordinates(L, -1);

  int rank = -1;
  lua_getfield(L, 1, "rank");
  if (lua_isinteger(L, -1))
    rank = (int)lua_tointeger(L, -1);

  polymesh_t* mesh = NULL;
  if (rank == -1) 
  {
    mesh = create_rectilinear_polymesh(comm, 
                                       xs->data, (int)xs->size,
                                       ys->data, (int)ys->size,
                                       zs->data, (int)zs->size);
  }
  else
  {
    mesh = create_rectilinear_polymesh_on_rank(comm, rank,
                                               xs->data, (int)xs->size,
                                               ys->data, (int)ys->size,
                                               zs->data, (int)zs->size);
  }

  lua_push_polymesh(L, mesh);
  return 1;
}

static lua_module_function polymeshes_funcs[] = {
  {"uniform", polymeshes_uniform, "polymeshes.uniform{comm = COMM, nx = NX, ny = NY, nz = NZ, bbox = BBOX} -> New uniform resolution mesh."},
  {"rectilinear", polymeshes_rectilinear, "polymeshes.rectilinear{comm = COMM, xs = XS, ys = YS, zs = ZS} -> New rectilinear mesh with defined node coordinates."},
  {NULL, NULL, NULL}
};

static lua_module_function points_funcs[] = {
  {NULL, NULL, NULL}
};

//------------------------------------------------------------------------
//                                API 
//------------------------------------------------------------------------

int lua_register_geometry_modules(lua_State* L)
{
  // Core types.
  lua_register_class(L, "coord_mapping", "A coordinate mapping.", cm_funcs, cm_methods);
  lua_register_class(L, "sd_func", "A signed distance function.", sd_funcs, sd_methods);
  lua_register_class(L, "sdt_func", "A time-dependent signed distance function.", sdt_funcs, sdt_methods);
  lua_register_record_type(L, "polymesh", "An arbitrary polyhedral mesh.", polymesh_funcs, polymesh_fields, polymesh_mm);

  // Register a module of mesh factory methods.
  lua_register_module(L, "polymeshes", "Functions for generating polymeshes.", polymeshes_funcs);

  // Register a module of point factory methods.
  lua_register_module(L, "points", "Functions for generating points.", points_funcs);

  return 0;
}

void lua_push_coord_mapping(lua_State* L, coord_mapping_t* X)
{
  lua_push_object(L, "coord_mapping", X, NULL);
}

bool lua_is_coord_mapping(lua_State* L, int index)
{
  return lua_is_object(L, index, "coord_mapping");
}

coord_mapping_t* lua_to_coord_mapping(lua_State* L, int index)
{
  return (coord_mapping_t*)lua_to_object(L, index, "coord_mapping");
}

void lua_push_sd_func(lua_State* L, sd_func_t* f)
{
  lua_push_object(L, "sd_func", f, NULL);
}

bool lua_is_sd_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sd_func");
}

sd_func_t* lua_to_sd_func(lua_State* L, int index)
{
  return (sd_func_t*)lua_to_object(L, index, "sd_func");
}

void lua_push_sdt_func(lua_State* L, sdt_func_t* f)
{
  lua_push_object(L, "sdt_func", f, NULL);
}

bool lua_is_sdt_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sdt_func");
}

sdt_func_t* lua_to_sdt_func(lua_State* L, int index)
{
  return (sdt_func_t*)lua_to_object(L, index, "sdt_func");
}

void lua_push_polymesh(lua_State* L, polymesh_t* m)
{
  lua_push_record(L, "polymesh", m, DTOR(polymesh_free));
}

bool lua_is_polymesh(lua_State* L, int index)
{
  return lua_is_record(L, index, "polymesh");
}

polymesh_t* lua_to_polymesh(lua_State* L, int index)
{
  return (polymesh_t*)lua_to_record(L, index, "polymesh");
}

