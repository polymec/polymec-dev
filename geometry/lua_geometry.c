// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "geometry/create_quad_planar_polymesh.h"
#include "geometry/create_hex_planar_polymesh.h"
#include "geometry/colmesh.h"

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

static void cm_map_vector(void* context, point_t* x, vector_t* v, vector_t* w)
{
  // Fetch our Lua table from the registry.
  lua_obj_t* lo = context;
  lua_pushlightuserdata(lo->L, lo);
  lua_gettable(lo->L, LUA_REGISTRYINDEX);

  // Call the map_vector method.
  lua_getfield(lo->L, -1, "map_vector");
  lua_pushvalue(lo->L, -2);
  lua_push_point(lo->L, x);
  lua_push_vector(lo->L, v);
  lua_call(lo->L, 3, 1);

  if (!lua_is_vector(lo->L, -1))
    luaL_error(lo->L, "map_vector method did not return a vector.");
  *w = *(lua_to_vector(lo->L, -1));
}

static void cm_compute_jacobian(void* context, point_t* x, real_t J[3][3])
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
  memcpy(J, lua_to_tensor2(lo->L, -1), sizeof(real_t) * 9);
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

  lua_getfield(L, 1, "map_point");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "map_point must be a method.");
  lua_setfield(L, obj_index, "map_point");
  vtable.map_point = cm_map_point;

  lua_getfield(L, 1, "map_vector");
  if (lua_isfunction(L, -1))
  {
    lua_setfield(L, obj_index, "map_vector");
    vtable.map_vector = cm_map_vector;
  }
  else if (!lua_isnil(L, -1))
    return luaL_error(L, "map_vector must be a method.");

  lua_getfield(L, 1, "jacobian");
  if (!lua_isfunction(L, -1))
    return luaL_error(L, "jacobian must be a method.");
  lua_setfield(L, obj_index, "jacobian");
  vtable.jacobian = cm_compute_jacobian;

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
  {"new", cm_new, "coord_mapping.new{name = NAME, map = Xp, [map_vector = Xv], jacobian = J} -> new coordinate mapping."},
  {"compose", cm_compose, "coord_mapping.compose(X1, X2) ­> Returns a coordinate mapping for X1 o X2."},
  {NULL, NULL, NULL}
};

static int cm_get_name(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  lua_pushstring(L, coord_mapping_name(X));
  return 1;
}

static int cm_get_inverse(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  coord_mapping_t* Xinv = coord_mapping_inverse(X);
  if (Xinv != NULL)
    lua_push_coord_mapping(L, Xinv);
  else
    lua_pushnil(L);
  return 1;
}

static int cm_set_inverse(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  coord_mapping_t* Xinv = lua_to_coord_mapping(L, 2);
  if (Xinv == NULL)
    luaL_error(L, "Argument must be a coord_mapping.");
  coord_mapping_set_inverse(X, Xinv);
  retain_ref(Xinv);
  return 0;
}

static lua_class_field cm_fields[] = {
  {"name", cm_get_name, NULL},
  {"inverse", cm_get_inverse, cm_set_inverse},
  {NULL, NULL, NULL}
};

static int cm_jacobian(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  if (X == NULL)
    luaL_error(L, "Method must be invoked with a coord_mapping.");
  if (!lua_is_point(L, 2))
    luaL_error(L, "Argument must be a point.");

  point_t* x = lua_to_point(L, 2);
  real_t J[3][3] = {{0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0},
                    {0.0, 0.0, 0.0}};
  coord_mapping_compute_jacobian(X, x, J);
  tensor2_t* J1 = tensor2_new(J[0][0], J[0][1], J[0][2],
                              J[1][0], J[1][1], J[1][2],
                              J[2][0], J[2][1], J[2][2]);
  lua_push_tensor2(L, J1);
  return 1;
}

static int cm_metric(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  if (X == NULL)
    luaL_error(L, "Method must be invoked with a coord_mapping.");
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

static int cm_call(lua_State* L)
{
  coord_mapping_t* X = lua_to_coord_mapping(L, 1);
  if (X == NULL)
    luaL_error(L, "Method must be invoked with a coord_mapping.");
  int num_args = lua_gettop(L);
  if (!(((num_args == 2) && lua_is_point(L, 2)) ||
        ((num_args == 3) && lua_is_point(L, 2) && lua_is_vector(L, 3))))
    luaL_error(L, "Argument must be a point or a point and a vector.");

  if (num_args == 2) // point
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
  if (X == NULL)
    luaL_error(L, "Method must be invoked with a coord_mapping.");
  lua_pushfstring(L, "coord_mapping '%s'", coord_mapping_name(X));
  return 1;
}

static lua_class_method cm_methods[] = {
  {"jacobian", cm_jacobian, "X:jacobian(x) -> Returns the jacobian of X at the point x."},
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

static int sd_get_name(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  lua_pushstring(L, sd_func_name(f));
  return 1;
}

static int sd_set_name(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sd_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static lua_class_field sd_fields[] = {
  {"name", sd_get_name, sd_set_name},
  {NULL, NULL, NULL}
};

static int sd_grad(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sd_func.");
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
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sd_func.");
  if (!lua_is_point(L, 2))
    return luaL_error(L, "Argument must be a point.");
  point_t* x = lua_to_point(L, 2);
  lua_pushnumber(L, sd_func_value(f, x));
  return 1;
}

static int sd_project(lua_State* L)
{
  sd_func_t* f = lua_to_sd_func(L, 1);
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sd_func.");
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
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sd_func.");
  lua_pushfstring(L, "sd_func '%s'", sd_func_name(f));
  return 1;
}

static lua_class_method sd_methods[] = {
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

static int sdt_get_name(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  lua_pushstring(L, sdt_func_name(f));
  return 1;
}

static int sdt_set_name(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sdt_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static lua_class_field sdt_fields[] = {
  {"name", sdt_get_name, sdt_set_name},
  {NULL, NULL, NULL}
};

static int sdt_grad(lua_State* L)
{
  sdt_func_t* f = lua_to_sdt_func(L, 1);
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sdt_func.");
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
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sdt_func.");
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
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sdt_func.");
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
  if (f == NULL)
    luaL_error(L, "Method must be invoked with an sdt_func.");
  lua_pushfstring(L, "sdt_func '%s'", sdt_func_name(f));
  return 1;
}

static lua_class_method sdt_methods[] = {
  {"grad", sdt_grad, "f:grad(x) -> Returns the gradient of f at x."},
  {"project", sdt_project, "f:project(x) -> Returns the projection of f to x."},
  {"__call", sdt_call, NULL},
  {"__tostring", sdt_tostring, NULL},
  {NULL, NULL, NULL}
};

static int p2_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args == 1)
  {
    if (!lua_istable(L, 1))
      return luaL_error(L, "Single argument must be a table with 3 or more vertices.");
  }
  else if (num_args < 3)
    return luaL_error(L, "Arguments must be vertices (at least 3).");

  point2_array_t* vertices = point2_array_new();
  if (num_args == 1)
  {
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, 1, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_is_point2(L, -1))
          luaL_error(L, "Item %d in argument is not a 2D point.", i);
        point2_array_append(vertices, *lua_to_point2(L, -1));
        lua_pop(L, 1);
      }
      ++i;
    }
  }
  else
  {
    for (int i = 1; i <= num_args; ++i)
    {
      if (!lua_is_point2(L, i))
        luaL_error(L, "Argument %d is not a 2D point.", i);
      point2_array_append(vertices, *lua_to_point2(L, i));
    }
  }
  if (vertices->size < 3)
    luaL_error(L, "A polygon must contain at least 3 vertices.");

  lua_push_polygon(L, polygon_new(vertices->data, vertices->size));
  point2_array_free(vertices);
  return 1;
}

static lua_module_function p2_funcs[] = {
  {"new", p2_new, "polygon.new(vertices) or polygon.new(v1, v2, ..., vn) -> new polygon in the plane with vertices given as 2D points."},
  {NULL, NULL, NULL}
};

static int p2_get_vertices(lua_State* L)
{
  polygon_t* p = lua_to_polygon(L, 1);
  lua_newtable(L);
  int pos = 0, i = 1;
  point2_t v;
  while (polygon_next_vertex(p, &pos, &v))
  {
    point2_t* v1 = point2_new(v.x, v.y);
    lua_push_point2(L, v1);
    lua_rawseti(L, -2, i);
    ++i;
  }
  return 1;
}

static int p2_get_area(lua_State* L)
{
  polygon_t* p = lua_to_polygon(L, 1);
  lua_push_real(L, polygon_area(p));
  return 1;
}

static int p2_get_centroid(lua_State* L)
{
  polygon_t* p = lua_to_polygon(L, 1);
  point2_t centroid;
  polygon_compute_centroid(p, &centroid);
  point2_t* xc = point2_new(centroid.x, centroid.y);
  lua_push_point2(L, xc);
  return 1;
}

static lua_class_field p2_fields[] = {
  {"vertices", p2_get_vertices, NULL},
  {"area", p2_get_area, NULL},
  {"centroid", p2_get_centroid, NULL},
  {NULL, NULL, NULL}
};

static int p2_tostring(lua_State* L)
{
  polygon_t* p = lua_to_polygon(L, 1);
  lua_pushfstring(L, "polygon (%d sides)", polygon_num_edges(p));
  if (p == NULL)
    luaL_error(L, "Method must be passed a polygon.");
  return 1;
}

static lua_class_method p2_methods[] = {
  {"__tostring", p2_tostring, NULL},
  {NULL, NULL, NULL}
};

static int p3_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_istable(L, 1) || !lua_istable(L, 2))
    return luaL_error(L, "Arguments are a list of vertices, and a table of faces consisting of lists of vertex indices.");

  // Get vertices.
  point_array_t* vertices = point_array_new();
  {
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, 1, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_is_point(L, -1))
          luaL_error(L, "Item %d in argument is not a point.", i);
        point_array_append(vertices, *lua_to_point(L, -1));
        lua_pop(L, 1);
      }
      ++i;
    }
  }
  if (vertices->size < 4)
    luaL_error(L, "A polyhedron must contain at least 4 vertices.");

  // Get faces.
  int_array_t* face_array = int_array_new();
  size_t_array_t* num_face_vertices = size_t_array_new();
  {
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, 2, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_istable(L, -1))
          luaL_error(L, "Item %d in argument 2 is not a list of vertex indices.", i);
        int j = 1;
        size_t num_vertices = 0;
        while (true)
        {
          lua_rawgeti(L, -1, j);
          if (lua_isnil(L, -1))
          {
            lua_pop(L, 1);
            break;
          }
          else
          {
            if (!lua_isinteger(L, -1))
              luaL_error(L, "Item %d in list of vertex indices for face %d is not an index.", j, i);
            int index = (int)lua_tointeger(L, -1);
            if ((index < 1) || (index > (int)(vertices->size)))
              luaL_error(L, "Vertex index %d for face %d is invalid: %d.", j, i, index);
            int_array_append(face_array, index);
            ++num_vertices;
            lua_pop(L, 1);
          }
          ++j;
        }
        if (num_vertices < 3)
          luaL_error(L, "Face %d must contain at least 3 vertices.", i);
        size_t_array_append(num_face_vertices, num_vertices);
        lua_pop(L, 1);
      }
      ++i;
    }
  }

  // Organize the face indices.
  int** faces = polymec_malloc(sizeof(int*) * num_face_vertices->size);
  size_t offset = 0;
  for (size_t f = 0; f < num_face_vertices->size; ++f)
  {
    size_t nv = num_face_vertices->data[f];
    faces[f] = polymec_malloc(sizeof(int) * nv);
    for (size_t fv = 0; fv < nv; ++fv, ++offset)
      faces[f][fv] = face_array->data[offset]-1;
  }
  int_array_free(face_array);

  // Create and push our polyhedron.
  lua_push_polyhedron(L, polyhedron_new(vertices->data, vertices->size,
                                        faces, num_face_vertices->data,
                                        num_face_vertices->size));

  // Clean up.
  for (size_t f = 0; f < num_face_vertices->size; ++f)
    polymec_free(faces[f]);
  polymec_free(faces);
  point_array_free(vertices);
  size_t_array_free(num_face_vertices);
  return 1;
}

static lua_module_function p3_funcs[] = {
  {"new", p3_new, "polyhedron.new(vertices, faces) -> new polyhedron with the given vertices and faces."},
  {NULL, NULL, NULL}
};

static int p3_get_vertices(lua_State* L)
{
  polyhedron_t* p = lua_to_polyhedron(L, 1);
  lua_newtable(L);
  int pos = 0, i = 1;
  point_t v;
  while (polyhedron_next_vertex(p, &pos, &v))
  {
    point_t* v1 = point_new(v.x, v.y, v.z);
    lua_push_point(L, v1);
    lua_rawseti(L, -2, i);
    ++i;
  }
  return 1;
}

static int p3_get_faces(lua_State* L)
{
  polyhedron_t* p = lua_to_polyhedron(L, 1);
  lua_newtable(L);
  int pos = 0, i = 1;
  point_t* face_vertices;
  size_t num_face_vertices;
  while (polyhedron_next_face(p, &pos, &face_vertices, &num_face_vertices, NULL, NULL))
  {
    lua_newtable(L);
    for (size_t v = 0; v < num_face_vertices; ++v)
    {
      point_t* v1 = point_new(face_vertices[v].x, face_vertices[v].y, face_vertices[v].z);
      lua_push_point(L, v1);
      lua_rawseti(L, -2, (int)(v+1));
    }
    lua_rawseti(L, -2, i);
    ++i;
  }
  return 1;
}

static int p3_get_volume(lua_State* L)
{
  polyhedron_t* p = lua_to_polyhedron(L, 1);
  lua_push_real(L, polyhedron_volume(p));
  return 1;
}

static int p3_get_centroid(lua_State* L)
{
  polyhedron_t* p = lua_to_polyhedron(L, 1);
  point_t centroid;
  polyhedron_compute_centroid(p, &centroid);
  point_t* xc = point_new(centroid.x, centroid.y, centroid.z);
  lua_push_point(L, xc);
  return 1;
}

static lua_class_field p3_fields[] = {
  {"vertices", p3_get_vertices, NULL},
  {"faces", p3_get_faces, NULL},
  {"volume", p3_get_volume, NULL},
  {"centroid", p3_get_centroid, NULL},
  {NULL, NULL, NULL}
};

static int p3_tostring(lua_State* L)
{
  polyhedron_t* p = lua_to_polyhedron(L, 1);
  lua_pushfstring(L, "polyhedron (%d faces)", polyhedron_num_faces(p));
  return 1;
}

static lua_class_method p3_methods[] = {
  {"__tostring", p3_tostring, NULL},
  {NULL, NULL, NULL}
};

static int pp_quad(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with nx, ny, nz, bbox fields.");

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

  lua_getfield(L, 1, "bbox");
  if (!lua_is_bbox(L, -1))
    luaL_error(L, "bbox must be a bounding box (bbox).");
  bbox_t* bbox = lua_to_bbox(L, -1);

  bool periodic_in_x = false;
  lua_getfield(L, 1, "periodic_in_x");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_x must be true or false.");
  periodic_in_x = lua_toboolean(L, -1);

  bool periodic_in_y = false;
  lua_getfield(L, 1, "periodic_in_y");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_y must be true or false.");
  periodic_in_y = lua_toboolean(L, -1);

  planar_polymesh_t* mesh = create_quad_planar_polymesh(nx, ny, bbox,
                                                        periodic_in_x,
                                                        periodic_in_y);
  lua_push_planar_polymesh(L, mesh);
  return 1;
}

static int pp_hex(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with radius and h fields.");

  lua_getfield(L, 1, "radius");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "radius must be a non-negative integer.");
  int radius = (int)lua_tointeger(L, -1);
  if (radius < 0)
    luaL_error(L, "radius must be non-negative.");

  lua_getfield(L, 1, "h");
  if (!lua_is_real(L, -1))
    luaL_error(L, "h must be a positive real number.");
  real_t h = lua_to_real(L, -1);
  if (h <= 0.0)
    luaL_error(L, "h must be positive.");

  planar_polymesh_t* mesh = create_hex_planar_polymesh(radius, h);
  lua_push_planar_polymesh(L, mesh);
  return 1;
}

static lua_module_function pp_funcs[] = {
  {"quad", pp_quad, "planar_polymesh.quad{nx = NX, ny = NY, bbox = BBOX, periodic_in_x = false, periodic_in_y = false} -> New uniform quadrilateral planar polymesh."},
  {"hex", pp_hex, "planar_polymesh.hex{nx = NX, ny = NY, bbox = BBOX, periodic_in_x = false, periodic_in_y = false} -> New uniform hexagonal planar polymesh."},
  {NULL, NULL, NULL}
};

static int pp_num_cells(lua_State* L)
{
  planar_polymesh_t* m = lua_to_planar_polymesh(L, 1);
  lua_newtable(L);
  lua_pushinteger(L, m->num_cells);
  return 1;
}

static int pp_num_edges(lua_State* L)
{
  planar_polymesh_t* m = lua_to_planar_polymesh(L, 1);
  lua_newtable(L);
  lua_pushinteger(L, m->num_edges);
  return 1;
}

static int pp_num_nodes(lua_State* L)
{
  planar_polymesh_t* m = lua_to_planar_polymesh(L, 1);
  lua_newtable(L);
  lua_pushinteger(L, m->num_nodes);
  return 1;
}

static lua_class_field pp_fields[] = {
  {"num_cells", pp_num_cells, NULL},
  {"num_edges", pp_num_edges, NULL},
  {"num_nodes", pp_num_nodes, NULL},
  {NULL, NULL, NULL}
};

static int pp_tostring(lua_State* L)
{
  planar_polymesh_t* m = lua_to_planar_polymesh(L, 1);
  lua_pushfstring(L, "planar_polymesh (%d cells)", m->num_cells);
  return 1;
}

static lua_class_method pp_methods[] = {
  {"__tostring", pp_tostring, NULL},
  {NULL, NULL, NULL}
};

static lua_module_function tagger_funcs[] = {
  {NULL, NULL, NULL}
};

static int t_create_tag(lua_State* L)
{
  tagger_t* t = lua_to_tagger(L, 1);
  if (t == NULL)
    return luaL_error(L, "Method must be invoked with a tagger.");
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument 1 must be a string.");
  if (!lua_istable(L, 3) && !lua_is_array(L, 3, LUA_ARRAY_INT))
    return luaL_error(L, "Argument 2 must be a table or array of integers.");
  int* tag = NULL;
  size_t size = 0;
  if (lua_is_array(L, 3, LUA_ARRAY_INT))
  {
    int_array_t* a = lua_to_array(L, 3, LUA_ARRAY_INT);
    tag = tagger_create_tag(t, lua_tostring(L, 2), a->size);
    size = a->size;
    memcpy(tag, a->data, sizeof(int) * size);
    int_array_free(a);
  }
  else
  {
    int_array_t* a = int_array_new();
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, 3, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_isinteger(L, -1))
          luaL_error(L, "Item %d in argument 1 is not an integer.", i);
        int_array_append(a, (int)(lua_tointeger(L, -1)));
        lua_pop(L, 1);
      }
      ++i;
    }
    tag = tagger_create_tag(t, lua_tostring(L, 2), a->size);
    size = a->size;
    memcpy(tag, a->data, sizeof(int) * size);
    int_array_free(a);
  }
  int_array_t* tt = int_array_new_with_data(tag, size);
  lua_push_array(L, tt, LUA_ARRAY_INT, true);
  return 1;
}

static int t_tag(lua_State* L)
{
  tagger_t* t = lua_to_tagger(L, 1);
  if (t == NULL)
    return luaL_error(L, "Method must be invoked with a tagger.");
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  size_t size;
  int* tag = tagger_tag(t, lua_tostring(L, 2), &size);
  if (tag != NULL)
  {
    int_array_t* tt = int_array_new_with_data(tag, size);
    lua_push_array(L, tt, LUA_ARRAY_INT, true);
    return 1;
  }
  else
    return 0;
}

static int t_has_tag(lua_State* L)
{
  tagger_t* t = lua_to_tagger(L, 1);
  if (t == NULL)
    return luaL_error(L, "Method must be invoked with a tagger.");
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  lua_pushboolean(L, tagger_has_tag(t, lua_tostring(L, 2)));
  return 1;
}

static int t_tostring(lua_State* L)
{
  lua_pushstring(L, "tagger");
  return 1;
}

static lua_class_method tagger_methods[] = {
  {"create_tag", t_create_tag, "tagger:create_tag(name, indices) -> Creates and returns a tag with the given name and indices."},
  {"tag", t_tag, "tagger:tag(name) -> Returns a tag with the given name."},
  {"has_tag", t_has_tag, "tagger:has_tag(name) -> Returns true if the tagger has a tag with the given name, false otherwise."},
  {"__tostring", t_tostring, NULL},
  {NULL, NULL, NULL}
};

static int pc_new(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) && (num_args != 3))
  {
    luaL_error(L, "Arguments must be an mpi.comm, a list of points, "
                  "and (optionally) a number of ghost points.");
  }

  if (!lua_is_mpi_comm(L, 1) )
    luaL_error(L, "Argument 1 must be an mpi.comm.");
  MPI_Comm comm = lua_to_mpi_comm(L, 1);

  point_t* points = NULL;
  int num_points = 0;
  bool free_points = false;
  if (!lua_istable(L, 2) && !lua_is_array(L, 2, LUA_ARRAY_POINT))
    luaL_error(L, "Argument 2 must be a list or an array of points.");
  if (lua_istable(L, 2))
  {
    point_array_t* points_array = point_array_new();
    int i = 1;
    while (true)
    {
      lua_rawgeti(L, 2, i);
      if (lua_isnil(L, -1))
      {
        lua_pop(L, 1);
        break;
      }
      else
      {
        if (!lua_is_point(L, -1))
          luaL_error(L, "Item %d in argument 1 is not a point.");
        point_t* x = lua_to_point(L, -1);
        point_array_append(points_array, *x);
        lua_pop(L, 1);
      }
      ++i;
    }
    points = points_array->data;
    num_points = (int)(points_array->size);
    point_array_release_data_and_free(points_array);
    free_points = true;
  }
  else
  {
    points = lua_to_array(L, 2, LUA_ARRAY_POINT);
    num_points = (int)lua_array_size(L, 2);
  }

  if ((num_args == 3) && !lua_isinteger(L, 3))
    luaL_error(L, "Argument 3 must be a number of ghost points.");
  int num_ghosts = (int)lua_tointeger(L, 3);
  if (num_ghosts < 0)
    luaL_error(L, "Number of ghost points must be positive.");

  // Create the point cloud.
  point_cloud_t* cloud = point_cloud_from_points(comm, points, (int)num_points);
  if (num_ghosts > 0)
    point_cloud_set_num_ghosts(cloud, num_ghosts);

  // Clean up.
  if (free_points)
    polymec_free(points);

  lua_push_point_cloud(L, cloud);
  return 1;
}

static int pc_repartition(lua_State* L)
{
  luaL_error(L, "can't repartition point clouds just yet!");
  return 0;
}

static lua_module_function pc_funcs[] = {
  {"new", pc_new, "point_cloud.new(comm, points [, num_ghosts]) -> new point cloud."},
  {"repartition", pc_repartition, "point_cloud.repartition(cloud) -> Repartitions the given point cloud."},
  {NULL, NULL, NULL}
};

static int pc_num_points(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushinteger(L, pc->num_points);
  return 1;
}

static int pc_num_ghosts(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushinteger(L, pc->num_ghosts);
  return 1;
}

static int pc_tags(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_push_tagger(L, pc->tags);
  return 1;
}

static lua_class_field pc_fields[] = {
  {"num_points", pc_num_points, NULL},
  {"num_ghosts", pc_num_ghosts, NULL},
  {"tags", pc_tags, NULL},
  {NULL, NULL, NULL}
};

static int pc_tostring(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushfstring(L, "point cloud (%d points)", pc->num_points);
  return 1;
}

static lua_class_method pc_methods[] = {
  {"__len", pc_num_points, NULL},
  {"__tostring", pc_tostring, NULL},
  {NULL, NULL, NULL}
};

static int um_new(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, bbox, npx, npy, npz, nx, ny, nz fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "bbox");
  if (!lua_is_bbox(L, -1))
    luaL_error(L, "bbox must be a bounding box.");
  bbox_t* bbox = lua_to_bbox(L, -1);
  if (bbox_is_empty_set(bbox))
    luaL_error(L, "bbox must be non-empty.");

  lua_getfield(L, 1, "npx");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "npx must be a positive integer.");
  int npx = (int)lua_tointeger(L, -1);
  if (npx < 1)
    luaL_error(L, "npx must be positive.");

  lua_getfield(L, 1, "npy");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "npy must be a positive integer.");
  int npy = (int)lua_tointeger(L, -1);
  if (npy < 1)
    luaL_error(L, "npy must be positive.");

  lua_getfield(L, 1, "npz");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "npz must be a positive integer.");
  int npz = (int)lua_tointeger(L, -1);
  if (npz < 1)
    luaL_error(L, "npz must be positive.");

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

  bool x_periodic = false;
  lua_getfield(L, 1, "periodic_in_x");
  if (lua_isboolean(L, -1))
    x_periodic = lua_toboolean(L, -1);
  else if (!lua_isnil(L, -1))
    luaL_error(L, "periodic_in_x must be a boolean.");

  bool y_periodic = false;
  lua_getfield(L, 1, "periodic_in_y");
  if (lua_isboolean(L, -1))
    y_periodic = lua_toboolean(L, -1);
  else if (!lua_isnil(L, -1))
    luaL_error(L, "periodic_in_y must be a boolean.");

  bool z_periodic = false;
  lua_getfield(L, 1, "periodic_in_z");
  if (lua_isboolean(L, -1))
    z_periodic = lua_toboolean(L, -1);
  else if (!lua_isnil(L, -1))
    luaL_error(L, "periodic_in_z must be a boolean.");

  unimesh_t* mesh = unimesh_new(comm, bbox, npx, npy, npz, nx, ny, nz,
                                x_periodic, y_periodic, z_periodic);
  lua_push_unimesh(L, mesh);
  return 1;
}

static int um_repartition(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_unimesh(L, 1))
    return luaL_error(L, "Argument must be a unimesh.");
  unimesh_t* mesh = lua_to_unimesh(L, 1);
  real_t imbalance_tol = 0.05;

  // Bug out if there's only one process.
  int nprocs;
  MPI_Comm_size(unimesh_comm(mesh), &nprocs);
  if (nprocs == 1)
    return 0;

  // Perform the repartitioning. FIXME: Add fields
  repartition_unimesh(&mesh, NULL, imbalance_tol, NULL, 0);
  return 0;
}

static lua_module_function um_funcs[] = {
  {"new", um_new, "unimesh.new{comm, bbox, npx, npy, npz, nx, ny, nz} -> "
                  "Creates a new uniform mesh on the communicator comm, which "
                  "fills the bounding box bbox with a lattice of npx x npy x npz "
                  "patches, each having nx x ny x nz identical cells."},
  {"repartition", um_repartition, "unimesh.repartition(m) -> Repartitions the colmesh m."},
  {NULL, NULL, NULL}
};

static int um_is_finalized(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  lua_pushboolean(L, unimesh_is_finalized(m));
  return 1;
}

static int um_bbox(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  lua_push_bbox(L, unimesh_bbox(m));
  return 1;
}

static int um_spacings(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  real_t dx, dy, dz;
  unimesh_get_spacings(m, &dx, &dy, &dz);
  lua_newtable(L);
  lua_push_real(L, dx);
  lua_rawseti(L, -2, 1);
  lua_push_real(L, dy);
  lua_rawseti(L, -2, 2);
  lua_push_real(L, dz);
  lua_rawseti(L, -2, 3);
  lua_pushvalue(L, -1);
  return 1;
}

static int um_extents(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  int npx, npy, npz;
  unimesh_get_extents(m, &npx, &npy, &npz);
  lua_newtable(L);
  lua_pushinteger(L, npx);
  lua_rawseti(L, -2, 1);
  lua_pushinteger(L, npy);
  lua_rawseti(L, -2, 2);
  lua_pushinteger(L, npz);
  lua_rawseti(L, -2, 3);
  lua_pushvalue(L, -1);
  return 1;
}

static int um_patch_size(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  int nx, ny, nz;
  unimesh_get_patch_size(m, &nx, &ny, &nz);
  lua_newtable(L);
  lua_pushinteger(L, nx);
  lua_rawseti(L, -2, 1);
  lua_pushinteger(L, ny);
  lua_rawseti(L, -2, 2);
  lua_pushinteger(L, nz);
  lua_rawseti(L, -2, 3);
  lua_pushvalue(L, -1);
  return 1;
}

static int um_patches(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  lua_newtable(L);
  int pos = 0, i, j, k, index = 1;
  while (unimesh_next_patch(m, &pos, &i, &j, &k, NULL))
  {
    lua_newtable(L);
    lua_pushinteger(L, i);
    lua_rawseti(L, -2, 1);
    lua_pushinteger(L, j);
    lua_rawseti(L, -2, 2);
    lua_pushinteger(L, k);
    lua_rawseti(L, -2, 3);
    lua_rawseti(L, -2, index);
    ++index;
  }
  lua_pushvalue(L, -1);
  return 1;
}

static lua_class_field um_fields[] = {
  {"is_finalized", um_is_finalized, NULL},
  {"bbox", um_bbox, NULL},
  {"spacings", um_spacings, NULL},
  {"extents", um_extents, NULL},
  {"patch_size", um_patch_size, NULL},
  {"patches", um_patches, NULL},
  {NULL, NULL, NULL}
};

static int um_tostring(lua_State* L)
{
  unimesh_t* m = lua_to_unimesh(L, 1);
  lua_pushfstring(L, "unimesh (%d local patches)",
                  unimesh_num_patches(m));
  return 1;
}

static lua_class_method um_methods[] = {
  {"__tostring", um_tostring, NULL},
  {NULL, NULL, NULL}
};

static int bm_new_(lua_State* L)
{
  if (!lua_is_mpi_comm(L, 1))
    luaL_error(L, "Argument 1 must be an mpi.comm.");
  MPI_Comm comm = lua_to_mpi_comm(L, 1);
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Argument 2 must be a positive number of x cells in a patch.");
  int patch_nx = (int)lua_tointeger(L, 2);
  if (patch_nx <= 0)
    luaL_error(L, "Argument 2 must be a positive number of x cells in a patch.");
  if (!lua_isinteger(L, 3))
    luaL_error(L, "Argument 3 must be a positive number of y cells in a patch.");
  int patch_ny = (int)lua_tointeger(L, 3);
  if (patch_ny <= 0)
    luaL_error(L, "Argument 2 must be a positive number of y cells in a patch.");
  if (!lua_isinteger(L, 4))
    luaL_error(L, "Argument 4 must be a positive number of z cells in a patch.");
  int patch_nz = (int)lua_tointeger(L, 4);
  if (patch_nz <= 0)
    luaL_error(L, "Argument 4 must be a positive number of z cells in a patch.");

  blockmesh_t* mesh = blockmesh_new(comm, patch_nx, patch_ny, patch_nz);
  lua_push_blockmesh(L, mesh);
  return 1;
}

static int bm_repartition(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_blockmesh(L, 1))
    return luaL_error(L, "Argument must be a blockmesh.");
  blockmesh_t* mesh = lua_to_blockmesh(L, 1);
  real_t imbalance_tol = 0.05;

  // Bug out if there's only one process.
  int nprocs;
  MPI_Comm_size(blockmesh_comm(mesh), &nprocs);
  if (nprocs == 1)
    return 0;

  // Perform the repartitioning. FIXME: Add fields
  repartition_blockmesh(&mesh, NULL, imbalance_tol, NULL, 0);
  return 0;
}

static lua_module_function bm_funcs[] = {
  {"new", bm_new_, "blockmesh.new(comm, patch_nx, patch_ny, patch_nz) -> New empty block-structured mesh with the given patch diemensions."},
  {"repartition", bm_repartition, "blockmesh.repartition(m) -> Repartitions the blockmesh m."},
  {NULL, NULL, NULL}
};

static int bm_is_finalized(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  lua_pushboolean(L, blockmesh_is_finalized(m));
  return 1;
}

static int bm_num_blocks(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  lua_pushinteger(L, (lua_Integer)blockmesh_num_blocks(m));
  return 1;
}

static int bm_comm(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  lua_push_mpi_comm(L, blockmesh_comm(m));
  return 1;
}

static lua_class_field bm_fields[] = {
  {"is_finalized", bm_is_finalized, NULL},
  {"num_blocks", bm_num_blocks, NULL},
  {"comm", bm_comm, NULL},
  {NULL, NULL, NULL}
};

static int bm_block(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  if (m == NULL)
    luaL_error(L, "Method must be invoked with a blockmesh.");
  if (!lua_isinteger(L, 2))
    return luaL_error(L, "Argument must be a block index.");
  int index = (int)(lua_tointeger(L, 2));
  if (index < 0)
    return luaL_error(L, "Block index must be non-negative.");
  if ((size_t)index >= blockmesh_num_blocks(m))
    return luaL_error(L, "Invalid block index: %d", index);
  lua_push_unimesh(L, blockmesh_block(m, index));
  return 1;
}

static int bm_add_block(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  if (m == NULL)
    luaL_error(L, "Method must be invoked with a blockmesh.");
  int num_args = lua_gettop(L);
  if (((num_args != 2) && (num_args == 4)) ||
      ((num_args == 2) && !lua_istable(L, 2)) ||
      ((num_args == 4) && (!lua_isinteger(L, 2) || !lua_isinteger(L, 3) ||
                           !lua_isinteger(L, 4))))
  {
    luaL_error(L, "mesh:add_block must be called with the following arguments, "
                  "or with a table with the following fields: "
                  "num_x_patches, num_y_patches, and num_z_patches");
  }
  bool table_args = (num_args == 2);

  int npx;
  if (table_args)
  {
    lua_getfield(L, 2, "num_x_patches");
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "num_x_patches must be a positive integer.");
    npx = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
  else
    npx = (int)lua_tointeger(L, 2);
  if (npx <= 0)
    return luaL_error(L, "num_x_patches must be positive.");

  int npy;
  if (table_args)
  {
    lua_getfield(L, 2, "num_y_patches");
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "num_y_patches must be a positive integer.");
    npy = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
  else
    npy = (int)lua_tointeger(L, 3);
  if (npy <= 0)
    return luaL_error(L, "num_y_patches must be positive.");

  int npz;
  if (table_args)
  {
    lua_getfield(L, 2, "num_z_patches");
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "num_z_patches must be a positive integer.");
    npz = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
  else
    npz = (int)lua_tointeger(L, 4);
  if (npz <= 0)
    return luaL_error(L, "num_z_patches must be positive.");

  blockmesh_add_block(m, npx, npy, npz);
  return 0;
}

static int bm_connect_blocks(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  if (m == NULL)
    luaL_error(L, "Method must be invoked with a blockmesh.");

  int num_args = lua_gettop(L);
  if (((num_args != 4) && (num_args != 2)) ||
      ((num_args == 2) && !lua_istable(L, 2)) ||
      ((num_args == 4) && (!lua_isinteger(L, 2) || !lua_istable(L, 3) ||
                           !lua_isinteger(L, 4) || !lua_istable(L, 5))))
  {
    luaL_error(L, "mesh:connect_blocks must be called with the following arguments, "
                  "or with a table containing these fields: "
                  "block1_index, block1_nodes, block2_index, and block2_nodes");
  }
  bool table_args = (num_args == 2);

  int index1;
  if (table_args)
  {
    lua_getfield(L, 2, "block1_index");
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "block1_index must be a valid block index.");
    index1 = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
  }
  else
    index1 = (int)lua_tointeger(L, 2);
  if ((index1 < 1) || ((size_t)index1 > blockmesh_num_blocks(m)))
    return luaL_error(L, "Invalid index for first block: %d.", index1);

  int nodes1[4], arg_index;
  if (table_args)
  {
    lua_getfield(L, 2, "block1_nodes");
    if (!lua_istable(L, -1))
      return luaL_error(L, "block1_nodes must be a list of 4 node indices.");
    arg_index = -1;
  }
  else
    arg_index = 3;
  for (int i = 1; i <= 4; ++i)
  {
    lua_rawgeti(L, arg_index, i);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "block1_nodes must be a list of 4 node indices.");
    int node = (int)lua_tointeger(L, -1);
    if ((node < 1) || (node > 8))
      return luaL_error(L, "block %d node %d is %d (must be between 1 and 8).", index1, i, node);
    nodes1[i-1] = node;
    lua_pop(L, 1);
  }
  if (table_args)
    lua_pop(L, 1);

  int index2;
  if (table_args)
  {
    lua_getfield(L, 2, "block2_index");
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "block2_index must be a valid block index.");
    index2 = (int)lua_tointeger(L, -1);
    if ((index2 < 1) || ((size_t)index2 > blockmesh_num_blocks(m)))
      return luaL_error(L, "Invalid index for second block: %d.", index2);
    lua_pop(L, 1);
  }
  else
    index2 = (int)lua_tointeger(L, 4);

  int nodes2[4];
  if (table_args)
  {
    lua_getfield(L, 2, "block2_nodes");
    if (!lua_istable(L, -1))
      return luaL_error(L, "block2_nodes must be a list of 4 node indices.");
    arg_index = -1;
  }
  else
    arg_index = 5;

  for (int i = 1; i <= 4; ++i)
  {
    lua_rawgeti(L, arg_index, i);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "block2_nodes must be a list of 4 node indices.");
    int node = (int)lua_tointeger(L, -1);
    if ((node < 1) || (node > 8))
      return luaL_error(L, "block %d node %d is %d (must be between 1 and 8).", index2, i, node);
    nodes2[i-1] = node;
    lua_pop(L, 1);
  }
  if (!table_args)
    lua_pop(L, 1);

  // Translate from Lua (1-based) indices to C (0-based) indices.
  int c_index1 = index1 - 1;
  int c_index2 = index2 - 1;
  int c_nodes1[4], c_nodes2[4];
  for (int i = 0; i < 4; ++i)
  {
    c_nodes1[i] = nodes1[i] - 1;
    c_nodes2[i] = nodes2[i] - 1;
  }

  char* err;
  if (!blockmesh_can_connect_blocks(m, c_index1, c_nodes1, c_index2, c_nodes2, &err))
    return luaL_error(L, "Couldn't connect blocks %d and %d: %s", index1, index2, err);

  blockmesh_connect_blocks(m, c_index1, c_nodes1, c_index2, c_nodes2);
  return 0;
}

static int bm_finalize(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  if (m == NULL)
    luaL_error(L, "Method must be invoked with a blockmesh.");
  if (blockmesh_is_finalized(m))
    luaL_error(L, "Blockmesh has already been finalized!");
  blockmesh_finalize(m);
  return 0;
}

static int bm_tostring(lua_State* L)
{
  blockmesh_t* m = lua_to_blockmesh(L, 1);
  if (m == NULL)
    luaL_error(L, "Method must be invoked with a blockmesh.");
  lua_pushfstring(L, "blockmesh (%d blocks)", (int)blockmesh_num_blocks(m));
  return 1;
}

static lua_class_method bm_methods[] = {
  {"block", bm_block, "mesh:block(index) -> Returns the block in the mesh "
                      "with the given index."},
  {"add_block", bm_add_block, "mesh:add_block(NPX, NPY, NPZ) OR "
                              "mesh:add_block{num_x_patches = NPX, "
                              "num_y_patches = NPY, num_z_patches = NPZ} "
                              "- Adds a new empty block to this mesh with the "
                              "given numbers of patches in x, y, and z."},
  {"connect_blocks", bm_connect_blocks, "mesh:connect_blocks(B1, B1_NODES, "
                                        "B2, B2_NODES) OR "
                                        "mesh:connect_blocks{block1_index = "
                                        "B1, block1_nodes = B1_NODES, "
                                        "block2_index = B2, "
                                        "block2_nodes = B2_NODES} -> Connects "
                                        "two blocks with indices B1 and B2 in "
                                        "the mesh, specifying the nodes to "
                                        "identify on the boundary with "
                                        "B1_NODES and B2_NODES."},
  {"finalize", bm_finalize, "mesh:finalize() - Finalizes a block mesh after "
                            "assembly."},
  {"__tostring", bm_tostring, NULL},
  {NULL, NULL, NULL}
};

static int colmesh_new_(lua_State* L)
{
  if (!lua_istable(L, 1))
  {
    luaL_error(L, "Argument must be a table with comm, columns, z1, z2, nz, "
                  "periodic_in_z fields.");
  }

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "columns");
  if (!lua_is_planar_polymesh(L, -1))
    luaL_error(L, "columns must be a planar polymesh.");
  planar_polymesh_t* columns = lua_to_planar_polymesh(L, -1);

  lua_getfield(L, 1, "nz");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "nz must be a positive integer.");
  int nz = (int)lua_tointeger(L, -1);
  if (nz < 1)
    luaL_error(L, "nz must be positive.");

  lua_getfield(L, 1, "z1");
  if (!lua_is_real(L, -1))
    luaL_error(L, "z1 must be a real number.");
  real_t z1 = lua_to_real(L, -1);

  lua_getfield(L, 1, "z2");
  if (!lua_is_real(L, -1))
    luaL_error(L, "z2 must be a real number.");
  real_t z2 = lua_to_real(L, -1);
  if (z2 <= z1)
    luaL_error(L, "z2 must be greater than z1.");

  bool periodic_in_z = false;
  lua_getfield(L, 1, "periodic_in_z");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_z must be true or false.");
  periodic_in_z = lua_toboolean(L, -1);

  colmesh_t* mesh = colmesh_new(comm, columns, z1, z2, nz, periodic_in_z);
  lua_push_colmesh(L, mesh);
  return 1;
}

static int colmesh_quad(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, nx, ny, nz, bbox fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm.");
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

  bool periodic_in_x = false;
  lua_getfield(L, 1, "periodic_in_x");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_x must be true or false.");
  periodic_in_x = lua_toboolean(L, -1);

  bool periodic_in_y = false;
  lua_getfield(L, 1, "periodic_in_y");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_y must be true or false.");
  periodic_in_y = lua_toboolean(L, -1);

  bool periodic_in_z = false;
  lua_getfield(L, 1, "periodic_in_z");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_z must be true or false.");
  periodic_in_z = lua_toboolean(L, -1);

  planar_polymesh_t* columns = create_quad_planar_polymesh(nx, ny, bbox,
                                                           periodic_in_x, periodic_in_y);
  colmesh_t* mesh = colmesh_new(comm, columns, bbox->z1, bbox->z2, nz, periodic_in_z);
  planar_polymesh_free(columns);
  lua_push_colmesh(L, mesh);
  return 1;
}

static int colmesh_hex(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, radius, h, z1, z2, nz fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm.");
  MPI_Comm comm = lua_to_mpi_comm(L, -1);

  lua_getfield(L, 1, "radius");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "radius must be a non-negative integer.");
  int radius = (int)lua_tointeger(L, -1);
  if (radius < 0)
    luaL_error(L, "radius must be non-negative.");

  lua_getfield(L, 1, "h");
  if (!lua_is_real(L, -1))
    luaL_error(L, "h must be a positive real number.");
  real_t h = lua_to_real(L, -1);
  if (h <= 0.0)
    luaL_error(L, "h must be positive.");

  lua_getfield(L, 1, "nz");
  if (!lua_isinteger(L, -1))
    luaL_error(L, "nz must be a positive integer.");
  int nz = (int)lua_tointeger(L, -1);
  if (nz < 1)
    luaL_error(L, "nz must be positive.");

  lua_getfield(L, 1, "z1");
  if (!lua_is_real(L, -1))
    luaL_error(L, "z1 must be a real number.");
  real_t z1 = lua_to_real(L, -1);

  lua_getfield(L, 1, "z2");
  if (!lua_is_real(L, -1))
    luaL_error(L, "z2 must be a real number.");
  real_t z2 = lua_to_real(L, -1);
  if (z2 <= z1)
    luaL_error(L, "z2 must be greater than z1.");

  bool periodic_in_z = false;
  lua_getfield(L, 1, "periodic_in_z");
  if (!lua_isnil(L, -1) && !lua_isboolean(L, -1))
    luaL_error(L, "periodic_in_z must be true or false.");
  periodic_in_z = lua_toboolean(L, -1);

  planar_polymesh_t* columns = create_hex_planar_polymesh(radius, h);
  colmesh_t* mesh = colmesh_new(comm, columns, z1, z2, nz, periodic_in_z);
  planar_polymesh_free(columns);
  lua_push_colmesh(L, mesh);
  return 1;
}

static int colmesh_repartition(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_colmesh(L, 1))
    return luaL_error(L, "Argument must be a colmesh.");
  colmesh_t* mesh = lua_to_colmesh(L, 1);
  real_t imbalance_tol = 0.05;

  // Bug out if there's only one process.
  int nprocs;
  MPI_Comm_size(colmesh_comm(mesh), &nprocs);
  if (nprocs == 1)
    return 0;

  // Perform the repartitioning. FIXME: Add fields
  repartition_colmesh(&mesh, NULL, imbalance_tol, NULL, 0);
  return 0;
}

static lua_module_function colmesh_funcs[] = {
  {"new", colmesh_new_, "colmesh.new{comm = COMM, columns = PLANAR_MESH, z1 = Z1, z2 = Z2, nz = NZ, periodic_in_z = false} -> New uniform quadrilateral prism mesh."},
  {"quad", colmesh_quad, "colmesh.quad{comm = COMM, bbox = BBOX, nx = NX, ny = NY, nz = NZ, periodic_in_x = false, periodic_in_y = false, periodic_in_z = false} -> New uniform quadrilateral prism mesh."},
  {"hex", colmesh_hex, "colmesh.hex{comm = COMM, radius = R, h = H, nz = NZ, z1 = Z1, z2 = Z2, periodic_in_z = false} -> New uniform hexagonal prism mesh."},
  {"repartition", colmesh_repartition, "colmesh.repartition(m) -> Repartitions the colmesh m."},
  {NULL, NULL, NULL}
};

static int cmesh_num_chunks(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  lua_pushinteger(L, colmesh_num_chunks(m));
  return 1;
}

static int cmesh_num_xy_chunks(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  colmesh_get_chunk_info(m, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  lua_pushinteger(L, (lua_Integer)num_xy_chunks);
  return 1;
}

static int cmesh_num_z_chunks(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  colmesh_get_chunk_info(m, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  lua_pushinteger(L, (lua_Integer)num_z_chunks);
  return 1;
}

static int cmesh_nz_per_chunk(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  colmesh_get_chunk_info(m, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  lua_pushinteger(L, (lua_Integer)nz_per_chunk);
  return 1;
}

static int cmesh_z1(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  real_t z1, z2;
  bool periodic;
  colmesh_get_z_info(m, &z1, &z2, &periodic);
  lua_push_real(L, z1);
  return 1;
}

static int cmesh_z2(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  real_t z1, z2;
  bool periodic;
  colmesh_get_z_info(m, &z1, &z2, &periodic);
  lua_push_real(L, z2);
  return 1;
}

static int cmesh_periodic_in_z(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  real_t z1, z2;
  bool periodic;
  colmesh_get_z_info(m, &z1, &z2, &periodic);
  lua_pushboolean(L, periodic);
  return 1;
}

static lua_class_field colmesh_fields[] = {
  {"num_chunks", cmesh_num_chunks, NULL},
  {"num_xy_chunks", cmesh_num_xy_chunks, NULL},
  {"num_z_chunks", cmesh_num_z_chunks, NULL},
  {"nz_per_chunk", cmesh_nz_per_chunk, NULL},
  {"z1", cmesh_z1, NULL},
  {"z2", cmesh_z2, NULL},
  {"periodic_in_z", cmesh_periodic_in_z, NULL},
  {NULL, NULL, NULL}
};

static int colmesh_tostring(lua_State* L)
{
  colmesh_t* m = lua_to_colmesh(L, 1);
  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  colmesh_get_chunk_info(m, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  real_t z1, z2;
  bool periodic;
  colmesh_get_z_info(m, &z1, &z2, &periodic);
  char periodic_str[12];
  if (periodic)
    snprintf(periodic_str, 11, " (periodic)");
  else
    periodic_str[0] = '\0';
  lua_pushfstring(L, "colmesh (%d chunks, %d xy chunks, %d z chunks, nz = %d, z1 = %g, z2 = %g%s)",
                  (int)colmesh_num_chunks(m), (int)num_xy_chunks, (int)num_z_chunks,
                  (int)nz_per_chunk, (double)z1, z2, periodic_str);
  return 1;
}

static lua_class_method colmesh_methods[] = {
  {"__tostring", colmesh_tostring, NULL},
  {NULL, NULL, NULL}
};

static int polymesh_uniform(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, nx, ny, nz, bbox fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm.");
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
  {
    mesh = create_uniform_polymesh(comm, nx, ny, nz, bbox);
    lua_push_polymesh(L, mesh);
  }
  else
  {
    mesh = create_uniform_polymesh_on_rank(comm, rank, nx, ny, nz, bbox);
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);
    if (rank == my_rank)
    {
      ASSERT(mesh != NULL);
      lua_push_polymesh(L, mesh);
    }
    else
    {
      ASSERT(mesh == NULL);
      lua_pushnil(L);
    }
  }
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

static int polymesh_rectilinear(lua_State* L)
{
  if (!lua_istable(L, 1))
    luaL_error(L, "Argument must be a table with comm, xs, ys, zs fields.");

  lua_getfield(L, 1, "comm");
  if (!lua_is_mpi_comm(L, -1))
    luaL_error(L, "comm must be an mpi.comm.");
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

  // Clean up.
  real_array_free(xs);
  real_array_free(ys);
  real_array_free(zs);

  lua_push_polymesh(L, mesh);
  return 1;
}

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

  // Perform the repartitioning. FIXME: Add fields
  bool result = repartition_polymesh(&mesh, NULL, imbalance_tol, NULL, 0);
  lua_pushboolean(L, result);
  return 1;
}

static lua_module_function polymesh_funcs[] = {
  {"uniform", polymesh_uniform, "polymesh.uniform{comm = COMM, nx = NX, ny = NY, nz = NZ, bbox = BBOX} -> New uniform resolution mesh."},
  {"rectilinear", polymesh_rectilinear, "polymesh.rectilinear{comm = COMM, xs = XS, ys = YS, zs = ZS} -> New rectilinear mesh with defined node coordinates."},
  {"repartition", polymesh_repartition, "mesh.repartition(m) -> Repartitions the polymesh m. Collective."},
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

static lua_class_field polymesh_fields[] = {
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

static lua_class_method polymesh_methods[] = {
  {"__tostring", polymesh_tostring, NULL},
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
  lua_register_class(L, "coord_mapping", "A coordinate mapping.", cm_funcs, cm_fields, cm_methods, NULL);
  lua_register_class(L, "sd_func", "A signed distance function.", sd_funcs, sd_fields, sd_methods, NULL);
  lua_register_class(L, "sdt_func", "A time-dependent signed distance function.", sdt_funcs, sdt_fields, sdt_methods, NULL);
  lua_register_class(L, "polygon", "A polygon in the plane.", p2_funcs, p2_fields, p2_methods, NULL);
  lua_register_class(L, "polyhedron", "A polyhedron.", p3_funcs, p3_fields, p3_methods, NULL);

  lua_register_class(L, "tagger", "An object that holds tags.", tagger_funcs, NULL, tagger_methods, NULL);
  lua_register_class(L, "point_cloud", "A point cloud in 3D space.", pc_funcs, pc_fields, pc_methods, DTOR(point_cloud_free));
  lua_register_class(L, "unimesh", "A uniform cartesian mesh.", um_funcs, um_fields, um_methods, DTOR(unimesh_free));
  lua_register_class(L, "blockmesh", "A block-structured mesh.", bm_funcs, bm_fields, bm_methods, DTOR(blockmesh_free));
  lua_register_class(L, "colmesh", "A polygonal extruded (pex) mesh.", colmesh_funcs, colmesh_fields, colmesh_methods, DTOR(colmesh_free));
  lua_register_class(L, "polymesh", "An arbitrary polyhedral mesh.", polymesh_funcs, polymesh_fields, polymesh_methods, DTOR(polymesh_free));
  lua_register_class(L, "planar_polymesh", "A planar polygonal mesh.", pp_funcs, pp_fields, pp_methods, DTOR(planar_polymesh_free));

  // Register some helpful symbols for unimesh.
  lua_getglobal(L, "unimesh");
  lua_pushinteger(L, (int)UNIMESH_X1_BOUNDARY);
  lua_setfield(L, -2, "x1_boundary");
  lua_pushinteger(L, (int)UNIMESH_X2_BOUNDARY);
  lua_setfield(L, -2, "x2_boundary");
  lua_pushinteger(L, (int)UNIMESH_Y1_BOUNDARY);
  lua_setfield(L, -2, "y1_boundary");
  lua_pushinteger(L, (int)UNIMESH_Y2_BOUNDARY);
  lua_setfield(L, -2, "y2_boundary");
  lua_pushinteger(L, (int)UNIMESH_Z1_BOUNDARY);
  lua_setfield(L, -2, "z1_boundary");
  lua_pushinteger(L, (int)UNIMESH_Z2_BOUNDARY);
  lua_setfield(L, -2, "z2_boundary");
  lua_pop(L, 1);

  // Register a module of point factory methods.
  lua_register_module(L, "points", "Functions for generating points.", NULL, points_funcs);

  return 0;
}

void lua_push_coord_mapping(lua_State* L, coord_mapping_t* X)
{
  lua_push_object(L, "coord_mapping", X);
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
  lua_push_object(L, "sd_func", f);
}

bool lua_is_sd_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sd_func");
}

sd_func_t* lua_to_sd_func(lua_State* L, int index)
{
  return (sd_func_t*)lua_to_object(L, index, "sd_func");
}

void lua_push_tagger(lua_State* L, tagger_t* t)
{
  lua_push_object(L, "tagger", t);
}

bool lua_is_tagger(lua_State* L, int index)
{
  return lua_is_object(L, index, "tagger");
}

tagger_t* lua_to_tagger(lua_State* L, int index)
{
  return (tagger_t*)lua_to_object(L, index, "tagger");
}

void lua_push_sdt_func(lua_State* L, sdt_func_t* f)
{
  lua_push_object(L, "sdt_func", f);
}

bool lua_is_sdt_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sdt_func");
}

sdt_func_t* lua_to_sdt_func(lua_State* L, int index)
{
  return (sdt_func_t*)lua_to_object(L, index, "sdt_func");
}

void lua_push_polygon(lua_State* L, polygon_t* p)
{
  lua_push_object(L, "polygon", p);
}

bool lua_is_polygon(lua_State* L, int index)
{
  return lua_is_object(L, index, "polygon");
}

polygon_t* lua_to_polygon(lua_State* L, int index)
{
  return (polygon_t*)lua_to_object(L, index, "polygon");
}

void lua_push_polyhedron(lua_State* L, polyhedron_t* p)
{
  lua_push_object(L, "polyhedron", p);
}

bool lua_is_polyhedron(lua_State* L, int index)
{
  return lua_is_object(L, index, "polyhedron");
}

polyhedron_t* lua_to_polyhedron(lua_State* L, int index)
{
  return (polyhedron_t*)lua_to_object(L, index, "polyhedron");
}

void lua_push_planar_polymesh(lua_State* L, planar_polymesh_t* m)
{
  lua_push_object(L, "planar_polymesh", m);
}

bool lua_is_planar_polymesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "planar_polymesh");
}

planar_polymesh_t* lua_to_planar_polymesh(lua_State* L, int index)
{
  return (planar_polymesh_t*)lua_to_object(L, index, "planar_polymesh");
}

void lua_push_point_cloud(lua_State* L, point_cloud_t* c)
{
  lua_push_object(L, "point_cloud", c);
}

bool lua_is_point_cloud(lua_State* L, int index)
{
  return lua_is_object(L, index, "point_cloud");
}

point_cloud_t* lua_to_point_cloud(lua_State* L, int index)
{
  return (point_cloud_t*)lua_to_object(L, index, "point_cloud");
}

void lua_push_unimesh(lua_State* L, unimesh_t* m)
{
  lua_push_object(L, "unimesh", m);
}

bool lua_is_unimesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "unimesh");
}

unimesh_t* lua_to_unimesh(lua_State* L, int index)
{
  return (unimesh_t*)lua_to_object(L, index, "unimesh");
}

void lua_push_blockmesh(lua_State* L, blockmesh_t* m)
{
  lua_push_object(L, "blockmesh", m);
}

bool lua_is_blockmesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "blockmesh");
}

blockmesh_t* lua_to_blockmesh(lua_State* L, int index)
{
  return (blockmesh_t*)lua_to_object(L, index, "blockmesh");
}

void lua_push_colmesh(lua_State* L, colmesh_t* m)
{
  lua_push_object(L, "colmesh", m);
}

bool lua_is_colmesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "colmesh");
}

colmesh_t* lua_to_colmesh(lua_State* L, int index)
{
  return (colmesh_t*)lua_to_object(L, index, "colmesh");
}

void lua_push_polymesh(lua_State* L, polymesh_t* m)
{
  lua_push_object(L, "polymesh", m);
}

bool lua_is_polymesh(lua_State* L, int index)
{
  return lua_is_object(L, index, "polymesh");
}

polymesh_t* lua_to_polymesh(lua_State* L, int index)
{
  return (polymesh_t*)lua_to_object(L, index, "polymesh");
}

