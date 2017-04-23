// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/lua_core.h"
#include "core/partition_mesh.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int z_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 2) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2))
  {
    return luaL_error(L, "Arguments must be real, imag components.");
  }

  real_t real = (real_t)lua_tonumber(L, 1);
  real_t imag = (real_t)lua_tonumber(L, 2);
  complex_t z = CMPLX(real, imag);
  lua_push_complex(L, z);
  return 1;
}

static int z_abs(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_pushnumber(L, (double)(cabs(z)));
  return 1;
}

static int z_arg(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_pushnumber(L, (double)(carg(z)));
  return 1;
}

static int z_conj(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_push_complex(L, conj(z));
  return 1;
}

static lua_module_function complex_funcs[] = {
  {"new", z_new},
  {"abs", z_abs},
  {"arg", z_arg},
  {"conj", z_conj},
  {NULL, NULL}
};

static int z_real(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_pushnumber(L, (double)creal(z));
  return 1;
}

static int z_imag(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_pushnumber(L, (double)cimag(z));
  return 1;
}

static lua_record_field complex_fields[] = {
  {"real", z_real, NULL},
  {"imag", z_imag, NULL},
  {NULL, NULL, NULL}
};

static int z_add(lua_State* L)
{
  complex_t z1 = lua_to_complex(L, 1);
  complex_t z2 = lua_to_complex(L, 2);
  lua_push_complex(L, z1 + z2);
  return 1;
}

static int z_sub(lua_State* L)
{
  complex_t z1 = lua_to_complex(L, 1);
  complex_t z2 = lua_to_complex(L, 2);
  lua_push_complex(L, z1 - z2);
  return 1;
}

static int z_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_complex(L, 2)) &&
      (!lua_is_complex(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a complex and a real.");
  complex_t z = lua_to_complex(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = (real_t)lua_tonumber(L, (lua_isnumber(L, 1)) ? 1 : 2);
  lua_push_complex(L, c * z);
  return 1;
}

static int z_div(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = (real_t)lua_tonumber(L, 2);
  lua_push_complex(L, z/c);
  return 1;
}

static int z_unm(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_push_complex(L, -z);
  return 1;
}

static int z_pow(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  complex_t p = 0.0;
  if (lua_is_complex(L, 2))
    p = lua_to_complex(L, 2);
  else if (lua_isnumber(L, 2))
    p = lua_tonumber(L, 2);
  else
    luaL_error(L, "Argument 2 must be a real or complex number.");
  lua_push_complex(L, cpow(z, p));
  return 1;
}

static int z_tostring(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  lua_pushfstring(L, "complex(%f, %f)", creal(z), cimag(z));
  return 1;
}

static lua_record_metamethod complex_mm[] = {
  {"__add", z_add},
  {"__sub", z_sub},
  {"__mul", z_mul},
  {"__div", z_div},
  {"__unm", z_unm},
  {"__pow", z_pow},
  {"__tostring", z_tostring},
  {NULL, NULL}
};

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

static int p_distance(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  point_t* q = lua_to_point(L, 2);
  if (q == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushnumber(L, (double)point_distance(p, q));
  return 1;
}

static lua_module_function point_funcs[] = {
  {"new", p_new},
  {"distance", p_distance},
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

static lua_record_field point_fields[] = {
  {"x", p_x, p_set_x},
  {"y", p_y, p_set_y},
  {"z", p_z, p_set_z},
  {NULL, NULL, NULL}
};

static int p_sub(lua_State* L)
{
  point_t* p1 = lua_to_point(L, 1);
  point_t* p2 = lua_to_point(L, 2);
  if (p2 == NULL)
    return luaL_error(L, "Arguments must both be points.");
  vector_t* diff = vector_new(p1->x - p2->x, p1->y - p2->y, p1->z - p2->z);
  lua_push_vector(L, diff);
  return 1;
}

static int p_len(lua_State* L)
{
  lua_pushnumber(L, 3.0);
  return 1;
}

static int p_tostring(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  lua_pushfstring(L, "point (%f, %f, %f)", p->x, p->y, p->z);
  return 1;
}

static lua_record_metamethod point_mm[] = {
  {"__sub", p_sub},
  {"__len", p_len},
  {"__tostring", p_tostring},
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

static int v_dot(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (w == NULL)
    return luaL_error(L, "Argument must be a vector.");
  lua_pushnumber(L, (double)vector_dot(v, w));
  return 1;
}

static lua_module_function vector_funcs[] = {
  {"new", v_new},
  {"dot", v_dot},
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

static lua_record_field vector_fields[] = {
  {"x", v_x, v_set_x},
  {"y", v_y, v_set_y},
  {"z", v_z, v_set_z},
  {NULL, NULL, NULL}
};

static int v_add(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (v == NULL)
    luaL_error(L, "Argument 1 must be a vector.");
  if (w == NULL)
    luaL_error(L, "Argument 2 must be a vector.");
  vector_t* sum = vector_new(v->x + w->x, v->y + w->y, v->z + w->z);
  lua_push_vector(L, sum);
  return 1;
}

static int v_sub(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (v == NULL)
    luaL_error(L, "Argument 1 must be a vector.");
  if (w == NULL)
    luaL_error(L, "Argument 2 must be a vector.");
  vector_t* diff = vector_new(v->x - w->x, v->y - w->y, v->z - w->z);
  lua_push_vector(L, diff);
  return 1;
}

static int v_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_vector(L, 2)) &&
      (!lua_is_vector(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a vector and a number.");
  vector_t* v = lua_to_vector(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = (real_t)lua_tonumber(L, (lua_isnumber(L, 1)) ? 1 : 2);
  vector_t* v1 = vector_new(c * v->x, c * v->y, c * v->z);
  lua_push_vector(L, v1);
  return 1;
}

static int v_div(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  if (v == NULL)
    luaL_error(L, "Argument 1 must be a vector.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = (real_t)lua_tonumber(L, 2);
  vector_t* v1 = vector_new(v->x/c, v->y/c, v->z/c);
  lua_push_vector(L, v1);
  return 1;
}

static int v_unm(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* v1 = vector_new(-1.0 * v->x, -1.0 * v->y, -1.0 * v->z);
  lua_push_vector(L, v1);
  return 1;
}

static int v_len(lua_State* L)
{
  lua_pushnumber(L, 3.0);
  return 1;
}

static int v_tostring(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushfstring(L, "vector (%f, %f, %f)", v->x, v->y, v->z);
  return 1;
}

static lua_record_metamethod vector_mm[] = {
  {"__add", v_add},
  {"__sub", v_sub},
  {"__mul", v_mul},
  {"__div", v_div},
  {"__unm", v_unm},
  {"__len", v_len},
  {"__tostring", v_tostring},
  {"dot", v_dot},
  {NULL, NULL}
};

static int t2_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 9) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2) || !lua_isnumber(L, 3) ||
      !lua_isnumber(L, 4) || !lua_isnumber(L, 5) || !lua_isnumber(L, 6) ||
      !lua_isnumber(L, 7) || !lua_isnumber(L, 8) || !lua_isnumber(L, 9))
  {
    return luaL_error(L, "Arguments must be xx, xy, xz, yx, yy, yz, zx, zy, zz components.");
  }

  real_t xx = (real_t)lua_tonumber(L, 1);
  real_t xy = (real_t)lua_tonumber(L, 2);
  real_t xz = (real_t)lua_tonumber(L, 3);
  real_t yx = (real_t)lua_tonumber(L, 4);
  real_t yy = (real_t)lua_tonumber(L, 5);
  real_t yz = (real_t)lua_tonumber(L, 6);
  real_t zx = (real_t)lua_tonumber(L, 7);
  real_t zy = (real_t)lua_tonumber(L, 8);
  real_t zz = (real_t)lua_tonumber(L, 9);
  tensor2_t* t = tensor2_new(xx, xy, xz,
                             yx, yy, yz,
                             zx, zy, zz);
  lua_push_tensor2(L, t);
  return 1;
}

static lua_module_function tensor2_funcs[] = {
  {"new", t2_new},
  {NULL, NULL}
};

static int t2_xx(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xx);
  return 1;
}

static int t2_set_xx(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xx = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_xy(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xy);
  return 1;
}

static int t2_set_xy(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xy = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_xz(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xz);
  return 1;
}

static int t2_set_xz(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_yx(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->yx);
  return 1;
}

static int t2_set_yx(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yx = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_yy(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->yy);
  return 1;
}

static int t2_set_yy(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yy = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_yz(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->yz);
  return 1;
}

static int t2_set_yz(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_zx(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->zx);
  return 1;
}

static int t2_set_zx(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->zx = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_zy(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->zy);
  return 1;
}

static int t2_set_zy(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->zy = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int t2_zz(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushnumber(L, (double)t->zz);
  return 1;
}

static int t2_set_zz(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->zz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static lua_record_field tensor2_fields[] = {
  {"xx", t2_xx, t2_set_xx},
  {"xy", t2_xy, t2_set_xy},
  {"xz", t2_xz, t2_set_xz},
  {"yx", t2_yx, t2_set_yx},
  {"yy", t2_yy, t2_set_yy},
  {"yz", t2_yz, t2_set_yz},
  {"zx", t2_zx, t2_set_zx},
  {"zy", t2_zy, t2_set_zy},
  {"zz", t2_zz, t2_set_zz},
  {NULL, NULL, NULL}
};

static int t2_add(lua_State* L)
{
  tensor2_t* t1 = lua_to_tensor2(L, 1);
  tensor2_t* t2 = lua_to_tensor2(L, 2);
  if (t1 == NULL)
    luaL_error(L, "Argument 1 must be a tensor2.");
  if (t2 == NULL)
    luaL_error(L, "Argument 2 must be a tensor2.");
  tensor2_t* sum = tensor2_new(t1->xx + t2->xx, t1->xy + t2->xy, t1->xz + t2->xz,
                               t1->yx + t2->yx, t1->yy + t2->yy, t1->yz + t2->yz,
                               t1->zx + t2->zx, t1->zy + t2->zy, t1->zz + t2->zz);
  lua_push_tensor2(L, sum);
  return 1;
}

static int t2_sub(lua_State* L)
{
  tensor2_t* t1 = lua_to_tensor2(L, 1);
  tensor2_t* t2 = lua_to_tensor2(L, 2);
  if (t1 == NULL)
    luaL_error(L, "Argument 1 must be a tensor2.");
  if (t2 == NULL)
    luaL_error(L, "Argument 2 must be a tensor2.");
  tensor2_t* diff = tensor2_new(t1->xx - t2->xx, t1->xy - t2->xy, t1->xz - t2->xz,
                                t1->yx - t2->yx, t1->yy - t2->yy, t1->yz - t2->yz,
                                t1->zx - t2->zx, t1->zy - t2->zy, t1->zz - t2->zz);
  lua_push_tensor2(L, diff);
  return 1;
}

static int t2_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_tensor2(L, 2)) &&
      (!lua_is_tensor2(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a tensor2 and a number.");
  tensor2_t* t = lua_to_tensor2(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = (real_t)lua_tonumber(L, (lua_isnumber(L, 1)) ? 1 : 2);
  tensor2_t* t1 = tensor2_new(c * t->xx, c * t->xy, c * t->xz,
                              c * t->yx, c * t->yy, c * t->yz,
                              c * t->zx, c * t->zy, c * t->zz);
  lua_push_tensor2(L, t1);
  return 1;
}

static int t2_div(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  if (t == NULL)
    luaL_error(L, "Argument 1 must be a tensor2.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = (real_t)lua_tonumber(L, 2);
  tensor2_t* t1 = tensor2_new(t->xx/c, t->xy/c, t->xz/c,
                              t->yx/c, t->yy/c, t->yz/c,
                              t->zx/c, t->zy/c, t->zz/c);
  lua_push_tensor2(L, t1);
  return 1;
}

static int t2_unm(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  tensor2_t* t1 = tensor2_new(-1.0 * t->xx, -1.0 * t->xy, -1.0 * t->xz,
                              -1.0 * t->yx, -1.0 * t->yy, -1.0 * t->yz,
                              -1.0 * t->zx, -1.0 * t->zy, -1.0 * t->zz);
  lua_push_tensor2(L, t1);
  return 1;
}

static int t2_len(lua_State* L)
{
  lua_pushnumber(L, 9.0);
  return 1;
}

static int t2_tostring(lua_State* L)
{
  tensor2_t* t = lua_to_tensor2(L, 1);
  lua_pushfstring(L, "tensor2 ((%f, %f, %f), (%f, %f, %f), (%f, %f, %f))", 
                  t->xx, t->xy, t->xz, t->yx, t->yy, t->yz, t->zx, t->zy, t->zz);
  return 1;
}

static lua_record_metamethod tensor2_mm[] = {
  {"__add", t2_add},
  {"__sub", t2_sub},
  {"__mul", t2_mul},
  {"__div", t2_div},
  {"__unm", t2_unm},
  {"__len", t2_len},
  {"__tostring", t2_tostring},
  {NULL, NULL}
};

static int st2_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 6) || 
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2) || !lua_isnumber(L, 3) ||
      !lua_isnumber(L, 4) || !lua_isnumber(L, 5) || !lua_isnumber(L, 6))
  {
    return luaL_error(L, "Arguments must be xx, xy, xz, yy, yz, zz components.");
  }

  real_t xx = (real_t)lua_tonumber(L, 1);
  real_t xy = (real_t)lua_tonumber(L, 2);
  real_t xz = (real_t)lua_tonumber(L, 3);
  real_t yy = (real_t)lua_tonumber(L, 4);
  real_t yz = (real_t)lua_tonumber(L, 5);
  real_t zz = (real_t)lua_tonumber(L, 6);
  sym_tensor2_t* t = sym_tensor2_new(xx, xy, xz,
                                         yy, yz,
                                             zz);
  lua_push_sym_tensor2(L, t);
  return 1;
}

static lua_module_function sym_tensor2_funcs[] = {
  {"new", st2_new},
  {NULL, NULL}
};

static int st2_xx(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xx);
  return 1;
}

static int st2_set_xx(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xx = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_xy(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xy);
  return 1;
}

static int st2_set_xy(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xy = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_xz(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xz);
  return 1;
}

static int st2_set_xz(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_yx(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xy);
  return 1;
}

static int st2_set_yx(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xy = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_yy(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->yy);
  return 1;
}

static int st2_set_yy(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yy = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_yz(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->yz);
  return 1;
}

static int st2_set_yz(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_zx(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->xz);
  return 1;
}

static int st2_set_zx(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_zy(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->yz);
  return 1;
}

static int st2_set_zy(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static int st2_zz(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushnumber(L, (double)t->zz);
  return 1;
}

static int st2_set_zz(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->zz = (real_t)lua_tonumber(L, 2);
  return 0;
}

static lua_record_field sym_tensor2_fields[] = {
  {"xx", st2_xx, st2_set_xx},
  {"xy", st2_xy, st2_set_xy},
  {"xz", st2_xz, st2_set_xz},
  {"yx", st2_yx, st2_set_yx},
  {"yy", st2_yy, st2_set_yy},
  {"yz", st2_yz, st2_set_yz},
  {"zx", st2_zx, st2_set_zx},
  {"zy", st2_zy, st2_set_zy},
  {"zz", st2_zz, st2_set_zz},
  {NULL, NULL, NULL}
};

static int st2_add(lua_State* L)
{
  sym_tensor2_t* t1 = lua_to_sym_tensor2(L, 1);
  sym_tensor2_t* t2 = lua_to_sym_tensor2(L, 2);
  if (t1 == NULL)
    luaL_error(L, "Argument 1 must be a tensor2.");
  if (t2 == NULL)
    luaL_error(L, "Argument 2 must be a tensor2.");
  sym_tensor2_t* sum = sym_tensor2_new(t1->xx + t2->xx, t1->xy + t2->xy, t1->xz + t2->xz,
                                                        t1->yy + t2->yy, t1->yz + t2->yz,
                                                                         t1->zz + t2->zz);
  lua_push_sym_tensor2(L, sum);
  return 1;
}

static int st2_sub(lua_State* L)
{
  sym_tensor2_t* t1 = lua_to_sym_tensor2(L, 1);
  sym_tensor2_t* t2 = lua_to_sym_tensor2(L, 2);
  if (t1 == NULL)
    luaL_error(L, "Argument 1 must be a tensor2.");
  if (t2 == NULL)
    luaL_error(L, "Argument 2 must be a tensor2.");
  sym_tensor2_t* diff = sym_tensor2_new(t1->xx - t2->xx, t1->xy - t2->xy, t1->xz - t2->xz,
                                                         t1->yy - t2->yy, t1->yz - t2->yz,
                                                         t1->zz - t2->zz);
  lua_push_sym_tensor2(L, diff);
  return 1;
}

static int st2_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_tensor2(L, 2)) &&
      (!lua_is_tensor2(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a tensor2 and a number.");
  sym_tensor2_t* t = lua_to_sym_tensor2(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = (real_t)lua_tonumber(L, (lua_isnumber(L, 1)) ? 1 : 2);
  sym_tensor2_t* t1 = sym_tensor2_new(c * t->xx, c * t->xy, c * t->xz,
                                                 c * t->yy, c * t->yz,
                                                            c * t->zz);
  lua_push_sym_tensor2(L, t1);
  return 1;
}

static int st2_div(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  if (t == NULL)
    luaL_error(L, "Argument 1 must be a tensor2.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = (real_t)lua_tonumber(L, 2);
  sym_tensor2_t* t1 = sym_tensor2_new(t->xx/c, t->xy/c, t->xz/c,
                                               t->yy/c, t->yz/c,
                                                        t->zz/c);
  lua_push_sym_tensor2(L, t1);
  return 1;
}

static int st2_unm(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  sym_tensor2_t* t1 = sym_tensor2_new(-1.0 * t->xx, -1.0 * t->xy, -1.0 * t->xz,
                                                    -1.0 * t->yy, -1.0 * t->yz,
                                                                  -1.0 * t->zz);
  lua_push_sym_tensor2(L, t1);
  return 1;
}

static int st2_len(lua_State* L)
{
  lua_pushnumber(L, 6.0);
  return 1;
}

static int st2_tostring(lua_State* L)
{
  sym_tensor2_t* t = lua_to_sym_tensor2(L, 1);
  lua_pushfstring(L, "sym_tensor2 ((%f, %f, %f), (%f, %f, %f), (%f, %f, %f))", 
                  t->xx, t->xy, t->xz, t->xy, t->yy, t->yz, t->xz, t->yz, t->zz);
  return 1;
}

static lua_record_metamethod sym_tensor2_mm[] = {
  {"__add", st2_add},
  {"__sub", st2_sub},
  {"__mul", st2_mul},
  {"__div", st2_div},
  {"__unm", st2_unm},
  {"__len", st2_len},
  {"__tostring", st2_tostring},
  {NULL, NULL}
};

static int bb_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a table containing x1, x2, y1, y2, z1, z2 values.");

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

static int bb_contains(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  point_t* p = lua_to_point(L, 2);
  if (p == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushboolean(L, bbox_contains(b, p));
  return 1;
}

static lua_module_function bbox_funcs[] = {
  {"new", bb_new},
  {NULL, NULL}
};

static int bb_tostring(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushfstring(L, "bbox (x1 = %f, x2 = %f, y1 = %f, y2 = %f, z1 = %f, z2 = %f)", 
                  b->x1, b->x2, b->y1, b->y2, b->z1, b->z2);
  return 1;
}

static lua_class_method bbox_methods[] = {
  {"contains", bb_contains},
  {"__tostring", bb_tostring},
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

static lua_module_function sp_funcs[] = {
  {"constant", sp_constant},
  {NULL, NULL}
};

static int sp_len(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  lua_pushinteger(L, sp_func_num_comp(f));
  return 1;
}

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

static lua_class_method sp_methods[] = {
  {"__len", sp_len},
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

static lua_module_function st_funcs[] = {
  {"constant", st_constant},
  {NULL, NULL}
};

static int st_len(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  lua_pushinteger(L, st_func_num_comp(f));
  return 1;
}

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

static lua_class_method st_methods[] = {
  {"__len", st_len},
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

static lua_module_function mesh_funcs[] = {
  {"repartition", mesh_repartition},
  {NULL, NULL}
};

static int mesh_num_cells(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_cells);
  return 1;
}

static int mesh_num_ghost_cells(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_ghost_cells);
  return 1;
}

static int mesh_num_faces(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_faces);
  return 1;
}

static int mesh_num_edges(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_edges);
  return 1;
}

static int mesh_num_nodes(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushinteger(L, m->num_nodes);
  return 1;
}

static lua_record_field mesh_fields[] = {
  {"num_cells", mesh_num_cells, NULL},
  {"num_ghost_cells", mesh_num_ghost_cells, NULL},
  {"num_faces", mesh_num_faces, NULL},
  {"num_edges", mesh_num_edges, NULL},
  {"num_nodes", mesh_num_nodes, NULL},
  {NULL, NULL, NULL}
};

static int mesh_tostring(lua_State* L)
{
  mesh_t* m = lua_to_mesh(L, 1);
  lua_pushfstring(L, "mesh (%d cells, %d faces, %d nodes)", 
                  m->num_cells, m->num_faces, m->num_nodes);
  return 1;
}

static lua_record_metamethod mesh_mm[] = {
  {"__tostring", mesh_tostring},
  {NULL, NULL}
};

static int pc_repartition(lua_State* L)
{
  luaL_error(L, "can't repartition point clouds just yet!");
  return 0;
}

static lua_module_function pc_funcs[] = {
  {"repartition", pc_repartition},
  {NULL, NULL}
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

static lua_record_field pc_fields[] = {
  {"num_points", pc_num_points, NULL},
  {"num_ghosts", pc_num_ghosts, NULL},
  {NULL, NULL, NULL}
};

static int pc_tostring(lua_State* L)
{
  point_cloud_t* pc = lua_to_point_cloud(L, 1);
  lua_pushfstring(L, "point cloud (%d points)", pc->num_points);
  return 1;
}

static lua_record_metamethod pc_mm[] = {
  {"__len", pc_num_points},
  {"__tostring", pc_tostring},
  {NULL, NULL}
};

static void lua_register_options(lua_State* L)
{
  // Create a new table and fill it with our named command line values.
  lua_newtable(L);
  options_t* opts = options_argv();
  int pos = 0;
  const char* opt_name;
  const char* opt_val;
  while (options_next_value(opts, &pos, &opt_name, &opt_val))
  {
    lua_pushstring(L, opt_val);
    lua_setfield(L, -2, opt_name);
  }
  lua_setglobal(L, "options");
}

//------------------------------------------------------------------------
//                                API 
//------------------------------------------------------------------------
extern int lua_register_array(lua_State* L); // defined in lua_array.c

int lua_register_core_modules(lua_State* L)
{
  // Core types.
  lua_register_record_type(L, "complex", complex_funcs, complex_fields, complex_mm);
  lua_register_record_type(L, "point", point_funcs, point_fields, point_mm);
  lua_register_record_type(L, "vector", vector_funcs, vector_fields, vector_mm);
  lua_register_record_type(L, "tensor2", tensor2_funcs, tensor2_fields, tensor2_mm);
  lua_register_record_type(L, "sym_tensor2", sym_tensor2_funcs, sym_tensor2_fields, sym_tensor2_mm);
  lua_register_array(L);

  lua_register_class(L, "bbox", bbox_funcs, bbox_methods);
  lua_register_class(L, "sp_func", sp_funcs, sp_methods);
  lua_register_class(L, "st_func", st_funcs, st_methods);
  lua_register_record_type(L, "mesh", mesh_funcs, mesh_fields, mesh_mm);
  lua_register_record_type(L, "point_cloud", pc_funcs, pc_fields, pc_mm);

  // Register the options table.
  lua_register_options(L);

  return 0;
}

void lua_push_complex(lua_State* L, complex_t z)
{
  complex_t* zz = polymec_malloc(sizeof(complex_t));
  *zz = z;
  lua_push_record(L, "complex", zz, polymec_free);
}

bool lua_is_complex(lua_State* L, int index)
{
  return lua_is_record(L, index, "complex");
}

complex_t lua_to_complex(lua_State* L, int index)
{
  return *((complex_t*)lua_to_record(L, index, "complex"));
}

void lua_push_point(lua_State* L, point_t* p)
{
  lua_push_record(L, "point", p, NULL);
}

bool lua_is_point(lua_State* L, int index)
{
  return lua_is_record(L, index, "point");
}

point_t* lua_to_point(lua_State* L, int index)
{
  return (point_t*)lua_to_record(L, index, "point");
}

void lua_push_vector(lua_State* L, vector_t* v)
{
  lua_push_record(L, "vector", v, NULL);
}

bool lua_is_vector(lua_State* L, int index)
{
  return lua_is_record(L, index, "vector");
}

vector_t* lua_to_vector(lua_State* L, int index)
{
  return (vector_t*)lua_to_record(L, index, "vector");
}

void lua_push_tensor2(lua_State* L, tensor2_t* t)
{
  lua_push_record(L, "tensor2", t, NULL);
}

bool lua_is_tensor2(lua_State* L, int index)
{
  return lua_is_record(L, index, "tensor2");
}

tensor2_t* lua_to_tensor2(lua_State* L, int index)
{
  return (tensor2_t*)lua_to_record(L, index, "tensor2");
}

void lua_push_sym_tensor2(lua_State* L, sym_tensor2_t* t)
{
  lua_push_record(L, "sym_tensor2", t, NULL);
}

bool lua_is_sym_tensor2(lua_State* L, int index)
{
  return lua_is_record(L, index, "sym_tensor2");
}

sym_tensor2_t* lua_to_sym_tensor2(lua_State* L, int index)
{
  return (sym_tensor2_t*)lua_to_record(L, index, "sym_tensor2");
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

