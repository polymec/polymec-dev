// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/lua_core.h"

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

  real_t real = lua_to_real(L, 1);
  real_t imag = lua_to_real(L, 2);
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
  {"new", z_new, "complex.new(re, im) -> new complex number."},
  {"abs", z_abs, "complex.abs(z) -> modulus of z."},
  {"arg", z_arg, "complex.arg(z) -> argument of z."},
  {"conj", z_conj, "complex.conj(z) -> complex conjugate of z."},
  {NULL, NULL, NULL}
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
  real_t c = lua_to_real(L, (lua_isnumber(L, 1)) ? 1 : 2);
  lua_push_complex(L, c * z);
  return 1;
}

static int z_div(lua_State* L)
{
  complex_t z = lua_to_complex(L, 1);
  if (!lua_is_real(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = lua_to_real(L, 2);
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
  else if (lua_is_real(L, 2))
    p = lua_to_real(L, 2);
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

  real_t x = lua_to_real(L, 1);
  real_t y = lua_to_real(L, 2);
  real_t z = lua_to_real(L, 3);
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
  {"new", p_new, "point.new(x, y, z) -> new 3D point."},
  {"distance", p_distance, "point.distance(x, y) -> |x - y|."},
  {NULL, NULL, NULL}
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
  p->x = lua_to_real(L, 2);
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
  p->y = lua_to_real(L, 2);
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
  p->z = lua_to_real(L, 2);
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
  lua_pushinteger(L, 3);
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

  real_t x = lua_to_real(L, 1);
  real_t y = lua_to_real(L, 2);
  real_t z = lua_to_real(L, 3);
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
  {"new", v_new, "vector.new(vx, vy, vz) -> new 3D vector."},
  {"dot", v_dot, "vector.dot(v1, v2) -> dot product of v1 with v2."},
  {NULL, NULL, NULL}
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
  v->x = lua_to_real(L, 2);
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
  v->y = lua_to_real(L, 2);
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
  v->z = lua_to_real(L, 2);
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
  real_t c = lua_to_real(L, (lua_isnumber(L, 1)) ? 1 : 2);
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
  real_t c = lua_to_real(L, 2);
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
  lua_pushinteger(L, 3);
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

  real_t xx = lua_to_real(L, 1);
  real_t xy = lua_to_real(L, 2);
  real_t xz = lua_to_real(L, 3);
  real_t yx = lua_to_real(L, 4);
  real_t yy = lua_to_real(L, 5);
  real_t yz = lua_to_real(L, 6);
  real_t zx = lua_to_real(L, 7);
  real_t zy = lua_to_real(L, 8);
  real_t zz = lua_to_real(L, 9);
  tensor2_t* t = tensor2_new(xx, xy, xz,
                             yx, yy, yz,
                             zx, zy, zz);
  lua_push_tensor2(L, t);
  return 1;
}

static lua_module_function tensor2_funcs[] = {
  {"new", t2_new, "tensor2.new(txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz) -> new general rank 2 tensor."},
  {NULL, NULL, NULL}
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
  t->xx = lua_to_real(L, 2);
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
  t->xy = lua_to_real(L, 2);
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
  t->xz = lua_to_real(L, 2);
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
  t->yx = lua_to_real(L, 2);
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
  t->yy = lua_to_real(L, 2);
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
  t->yz = lua_to_real(L, 2);
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
  t->zx = lua_to_real(L, 2);
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
  t->zy = lua_to_real(L, 2);
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
  t->zz = lua_to_real(L, 2);
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
  real_t c = lua_to_real(L, (lua_isnumber(L, 1)) ? 1 : 2);
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
  real_t c = lua_to_real(L, 2);
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
  lua_pushinteger(L, 9);
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

  real_t xx = lua_to_real(L, 1);
  real_t xy = lua_to_real(L, 2);
  real_t xz = lua_to_real(L, 3);
  real_t yy = lua_to_real(L, 4);
  real_t yz = lua_to_real(L, 5);
  real_t zz = lua_to_real(L, 6);
  sym_tensor2_t* t = sym_tensor2_new(xx, xy, xz,
                                         yy, yz,
                                             zz);
  lua_push_sym_tensor2(L, t);
  return 1;
}

static lua_module_function sym_tensor2_funcs[] = {
  {"new", st2_new, "sym_tensor2.new(txx, txy, txz, tyy, tyz, tzz) -> new symmetric rank 2 tensor."},
  {NULL, NULL, NULL}
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
  t->xx = lua_to_real(L, 2);
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
  t->xy = lua_to_real(L, 2);
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
  t->xz = lua_to_real(L, 2);
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
  t->xy = lua_to_real(L, 2);
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
  t->yy = lua_to_real(L, 2);
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
  t->yz = lua_to_real(L, 2);
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
  t->xz = lua_to_real(L, 2);
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
  t->yz = lua_to_real(L, 2);
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
  t->zz = lua_to_real(L, 2);
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
    luaL_error(L, "Argument 1 must be a sym_tensor2.");
  if (t2 == NULL)
    luaL_error(L, "Argument 2 must be a sym_tensor2.");
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
    luaL_error(L, "Argument 1 must be a sym_tensor2.");
  if (t2 == NULL)
    luaL_error(L, "Argument 2 must be a sym_tensor2.");
  sym_tensor2_t* diff = sym_tensor2_new(t1->xx - t2->xx, t1->xy - t2->xy, t1->xz - t2->xz,
                                                         t1->yy - t2->yy, t1->yz - t2->yz,
                                                         t1->zz - t2->zz);
  lua_push_sym_tensor2(L, diff);
  return 1;
}

static int st2_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_sym_tensor2(L, 2)) &&
      (!lua_is_sym_tensor2(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a sym_tensor2 and a number.");
  sym_tensor2_t* t = lua_to_sym_tensor2(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = lua_to_real(L, (lua_isnumber(L, 1)) ? 1 : 2);
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
    luaL_error(L, "Argument 1 must be a sym_tensor2.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = lua_to_real(L, 2);
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
  lua_pushinteger(L, 6);
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

static int mpi_comm_split(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    return luaL_error(L, "Arguments must be an mpi.comm, a color, and a key.");
  if (!lua_is_mpi_comm(L, 1))
    return luaL_error(L, "Argument 1 must be an mpi.comm.");
  if (!lua_isinteger(L, 2))
    return luaL_error(L, "Argument 2 must be an integer color.");
  if (!lua_isinteger(L, 3))
    return luaL_error(L, "Argument 3 must be an integer key.");
  MPI_Comm comm = lua_to_mpi_comm(L, 1);
  int color = (int)lua_tointeger(L, 2);
  int key = (int)lua_tointeger(L, 3);

  MPI_Comm new_comm;
  MPI_Comm_split(comm, color, key, &new_comm);
  lua_push_mpi_comm(L, new_comm);

  return 1;
}

static lua_module_function mpi_comm_funcs[] = {
  {"split", mpi_comm_split, "mpi.comm.split(comm, color, key) -> new communicator."},
  {NULL, NULL, NULL}
};

static int mpi_comm_size(lua_State* L)
{
  MPI_Comm comm = lua_to_mpi_comm(L, 1);
  int size;
  MPI_Comm_size(comm, &size);
  lua_pushinteger(L, size);
  return 1;
}

static int mpi_comm_rank(lua_State* L)
{
  MPI_Comm comm = lua_to_mpi_comm(L, 1);
  int rank;
  MPI_Comm_rank(comm, &rank);
  lua_pushinteger(L, rank);
  return 1;
}

static lua_record_field mpi_comm_fields[] = {
  {"size", mpi_comm_size, NULL},
  {"rank", mpi_comm_rank, NULL},
  {NULL, NULL, NULL}
};

static int mpi_comm_tostring(lua_State* L)
{
  MPI_Comm comm = lua_to_mpi_comm(L, 1);
  if (comm == MPI_COMM_WORLD)
    lua_pushstring(L, "mpi.comm (MPI_COMM_WORLD)");
  else if (comm == MPI_COMM_SELF)
    lua_pushstring(L, "mpi.comm (MPI_COMM_SELF)");
  else
  {
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    lua_pushfstring(L, "mpi.comm (%d procs)", nprocs);
  }
  return 1;
}

static lua_record_metamethod mpi_comm_mm[] = {
  {"__tostring", mpi_comm_tostring},
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
    real_t num = (i % 2) ? -REAL_MAX : REAL_MAX;
    if (lua_isnumber(L, -1))
      num = lua_to_real(L, -1);
    else if (!lua_isnil(L, -1) && !lua_isnumber(L, -1))
    {
      return luaL_error(L, "Invalid entry for '%s'.\n"
                        "x1, x2, y1, y2, z1, z2, must all be numbers.", entries[i]);
    }
    switch(i)
    {
      case 0: bbox->x1 = num;
              break;
      case 1: bbox->x2 = num;
              break;
      case 2: bbox->y1 = num;
              break;
      case 3: bbox->y2 = num;
              break;
      case 4: bbox->z1 = num;
              break;
      case 5: bbox->z2 = num;
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
  {"new", bb_new, "bbox.new{x1 = X1, x2 = X2, y1 = Y1, y2 = Y2, z1 = Z1, z2 = Z2} -> new bounding box."},
  {NULL, NULL, NULL}
};

static int bb_tostring(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  lua_pushfstring(L, "bbox (x1 = %f, x2 = %f, y1 = %f, y2 = %f, z1 = %f, z2 = %f)", 
                  b->x1, b->x2, b->y1, b->y2, b->z1, b->z2);
  return 1;
}

static lua_class_method bbox_methods[] = {
  {"contains", bb_contains, "box.contains(x) -> true if box contains x, false otherwise."},
  {"__tostring", bb_tostring, NULL},
  {NULL, NULL, NULL}
};

static int sp_constant(lua_State* L)
{
  // Check the argument.
  int num_args = lua_gettop(L);
  real_t val[num_args];
  for (int i = 1; i <= num_args; ++i)
  {
    if (!lua_isnumber(L, i))
      return luaL_error(L, "Argument %d must be a number.", i);
    val[i-1] = lua_to_real(L, i);
  }
  sp_func_t* f = constant_sp_func_new(val, num_args);
  lua_push_sp_func(L, f);
  return 1;
}

static lua_module_function sp_funcs[] = {
  {"constant", sp_constant, "sp_func.constant(F0) -> returns a function with the constant value F0."},
  {NULL, NULL, NULL}
};

static int sp_rename(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sp_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static int sp_gc(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  f = NULL;
  return 0;
}

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

static int sp_tostring(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  char homo_str[33];
  if (sp_func_is_homogeneous(f))
    strcpy(homo_str, "homogeneous, ");
  else
    homo_str[0] = '\0';
  lua_pushfstring(L, "sp_func '%s' (%s%d components)", 
                  sp_func_name(f), homo_str, sp_func_num_comp(f));
  return 1;
}

static lua_class_method sp_methods[] = {
  {"rename", sp_rename, "f:rename(name) -> renames f to the given name."},
  {"__gc", sp_gc, NULL},
  {"__len", sp_len, NULL},
  {"__call", sp_call, NULL},
  {"__tostring", sp_tostring, NULL},
  {NULL, NULL, NULL}
};

static int st_constant(lua_State* L)
{
  // Check the argument.
  int num_args = lua_gettop(L);
  real_t val[num_args];
  for (int i = 1; i <= num_args; ++i)
  {
    if (!lua_isnumber(L, i))
      return luaL_error(L, "Argument %d must be a number.", i);
    val[i-1] = lua_to_real(L, i);
  }
  st_func_t* f = constant_st_func_new(val, num_args);
  lua_push_st_func(L, f);
  return 1;
}

static int st_from_sp_func(lua_State* L)
{
  // Check the argument.
  if (!lua_is_sp_func(L, 1))
    return luaL_error(L, "Argument must be an sp_func.");
  sp_func_t* f = lua_to_sp_func(L, 1);
  st_func_t* g = st_func_from_sp_func(f);
  lua_push_st_func(L, g);
  return 1;
}

static lua_module_function st_funcs[] = {
  {"constant", st_constant, "st_func.constant(F0) -> returns a function with the constant value F0."},
  {"from_sp_func", st_from_sp_func, "st_func.from_sp_func(f) -> returns a time-dependent function identical to f."},
  {NULL, NULL, NULL}
};

static int st_rename(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  st_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static int st_freeze(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Argument must be a time.");
  real_t t = lua_to_real(L, 3);
  sp_func_t* g = st_func_freeze(f, t);
  lua_push_sp_func(L, g);
  return 1;
}

static int st_gc(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  f = NULL;
  return 0;
}

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
  if (!lua_isnumber(L, 3))
    return luaL_error(L, "Second argument must be a time.");
  point_t* x = lua_to_point(L, 2);
  real_t t = lua_to_real(L, 3);
  int nc = st_func_num_comp(f);
  real_t val[nc];
  st_func_eval(f, x, t, val);
  for (int i = 0; i < nc; ++i)
    lua_pushnumber(L, val[i]);
  return nc;
}

static int st_tostring(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  char homo_str[33];
  if (st_func_is_homogeneous(f))
    strcpy(homo_str, "homogeneous, ");
  else
    homo_str[0] = '\0';
  char const_str[33];
  if (st_func_is_constant(f))
    strcpy(const_str, "constant, ");
  else
    const_str[0] = '\0';
  lua_pushfstring(L, "sp_func '%s' (%s%s%d components)", 
                  st_func_name(f), homo_str, const_str, st_func_num_comp(f));
  return 1;
}

static lua_class_method st_methods[] = {
  {"rename", st_rename, "f:rename(name) -> renames f to the given name."},
  {"freeze", st_freeze, "f:freeze(t) -> freezes value of f at time t."},
  {"__gc", st_gc, NULL},
  {"__len", st_len, NULL},
  {"__call", st_call, NULL},
  {"__tostring", st_tostring, NULL},
  {NULL, NULL, NULL}
};

static lua_module_function tagger_funcs[] = {
  {NULL, NULL, NULL}
};

static int t_create_tag(lua_State* L)
{
  tagger_t* t = lua_to_tagger(L, 1);
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
  }
  else
  {
    int_array_t* a = int_array_new();
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
        if (!lua_isinteger(L, -1))
          luaL_error(L, "Item %d in argument 1 is not an integer.");
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
    luaL_error(L, "Arguments must be an MPI communicator, a list of points, "
                  "and (optionally) a number of ghost points.");
  }

  if (!lua_is_mpi_comm(L, 1) )
    luaL_error(L, "Argument 1 must be an MPI communicator.");
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

static lua_record_field pc_fields[] = {
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

// Here's a dir() function just like Python's.
static int lua_dir(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    luaL_error(L, "dir() accepts exactly one argument.");

  int table_index = 0;
  if (lua_istable(L, 1))
    table_index = 1;
  else if (lua_isuserdata(L, 1) && (lua_getmetatable(L, 1) != 0))
    table_index = 2;
  else
  {
    // Return an empty table.
    lua_newtable(L);
    return 1;
  }

  // Go through the keys of the table and push them into a new table.
  lua_newtable(L);
  int index = 0;
  lua_pushnil(L);
  while (lua_next(L, table_index))
  {
    ++index;

    // Key is at index -2, value is at -1, new table is at -3.
    // First, replace the value with a copy of the key.
    lua_pop(L, 1);
    lua_pushvalue(L, -1);

    // Now add this key to our new table at -3. This pops the key copy 
    // off the stack.
    lua_rawseti(L, -3, index); 
  }

  // Clean up the stack.
  if (table_index == 2)
    lua_remove(L, table_index);

  ASSERT(lua_gettop(L) == num_args + 1);
  return 1;
}

static int lua_do_nothing(lua_State* L)
{
  // This is just a stub.
  return 0;
}

static int fake_io_open(lua_State* L)
{
  // We just return a table with some benign entries.
  lua_newtable(L);
  lua_pushcfunction(L, lua_do_nothing);
  lua_setfield(L, -2, "read");
  lua_pushcfunction(L, lua_do_nothing);
  lua_setfield(L, -2, "write");
  lua_pushcfunction(L, lua_do_nothing);
  lua_setfield(L, -2, "close");
  return 1;
}

extern void lua_get_docstrings(lua_State* L);
static int lua_help(lua_State* L)
{
  // Dig up our docstrings table.
  lua_get_docstrings(L);

  // Retrieve the record for the argument.
  lua_pushnil(L);
  lua_copy(L, -3, -1);
  lua_gettable(L, -2);
  if (lua_isstring(L, -1))
  {
    const char* doc = lua_tostring(L, -1);
    fprintf(stdout, "%s\n", doc);
  }
  lua_pop(L, 2);
  return 0;
}

extern void lua_set_docstring(lua_State* L, int index, const char* docstring);
static int lua_document(lua_State* L)
{
  if (lua_isnil(L, 1))
    luaL_error(L, "Argument 1 cannot be nil.");
  if (lua_isstring(L, 2))
    lua_set_docstring(L, 1, lua_tostring(L, 2));
  else
    luaL_error(L, "Argument 2 must be a docstring for argument 1.");
  return 0;
}

extern void lua_replace_tostring(lua_State* L);

static void lua_register_util_funcs(lua_State* L)
{
  // Python-like dir() function.
  lua_pushcfunction(L, lua_dir);
  lua_setglobal(L, "dir");
  lua_getglobal(L, "dir");
  lua_set_docstring(L, -1, "dir(X) -> Shows the contents of X (if it's a table or a module).");
  lua_pop(L, 1);

  // Python-like help() function.
  lua_pushcfunction(L, lua_help);
  lua_setglobal(L, "help");
  lua_getglobal(L, "help");
  lua_set_docstring(L, -1, "help(X) -> Prints help on X, which can a module, a function, or an object.");
  lua_pop(L, 1);

  // A string for setting a docstring.
  lua_pushcfunction(L, lua_document);
  lua_setglobal(L, "document");
  lua_getglobal(L, "document");
  lua_set_docstring(L, -1, "document(X, S) -> Documents the object X with the docstring S.");
  lua_pop(L, 1);

  lua_getglobal(L, "print");
  lua_set_docstring(L, -1, "print(X, ...) -> Prints the value of X (etc) to the screen.");
  lua_pop(L, 1);

  // Replace the default tostring function.
  lua_replace_tostring(L);

  // On MPI ranks != 0, neutralize I/O functions.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0)
  {
    // print().
    lua_pushcfunction(L, lua_do_nothing);
    lua_setglobal(L, "print");

    // help().
    lua_pushcfunction(L, lua_do_nothing);
    lua_setglobal(L, "help");

    // io module.
    lua_getglobal(L, "io");
    lua_pushcfunction(L, fake_io_open);
    lua_setfield(L, -2, "open");
    lua_pushcfunction(L, lua_do_nothing);
    lua_setfield(L, -2, "close");
    lua_pop(L, 1);
  }
}

//------------------------------------------------------------------------
//                                API 
//------------------------------------------------------------------------
extern int lua_register_array(lua_State* L); // defined in lua_array.c
extern int lua_register_ndarray(lua_State* L); // defined in lua_ndarray.c

static int lua_register_constants(lua_State* L)
{
  // Register the origin.
  lua_getglobal(L, "point");
  point_t* O = point_new(0.0, 0.0, 0.0);
  lua_push_point(L, O);
  lua_setfield(L, -2, "origin");
  O = NULL;
  lua_pop(L, 1);

  // Register unit vectors in 3D cartesian coordinates.
  lua_getglobal(L, "vector");
  vector_t* e1 = vector_new(1.0, 0.0, 0.0);
  lua_push_vector(L, e1);
  lua_setfield(L, -2, "e1");
  e1 = NULL;
  vector_t* e2 = vector_new(0.0, 1.0, 0.0);
  lua_push_vector(L, e2);
  lua_setfield(L, -2, "e2");
  e2 = NULL;
  vector_t* e3 = vector_new(0.0, 0.0, 1.0);
  lua_push_vector(L, e3);
  lua_setfield(L, -2, "e3");
  e3 = NULL;
  lua_pop(L, 1);

  // Register unit tensors.
  lua_getglobal(L, "tensor2");
  tensor2_t I = {.xx = 1.0, .xy = 0.0, .xz = 0.0,
                 .yx = 0.0, .yy = 1.0, .yz = 0.0,
                 .zx = 0.0, .zy = 0.0, .zz = 1.0};
  lua_push_tensor2(L, &I);
  lua_setfield(L, -2, "unit");
  lua_pop(L, 1);

  lua_getglobal(L, "sym_tensor2");
  sym_tensor2_t I_sym = {.xx = 1.0, .xy = 0.0, .xz = 0.0,
                                    .yy = 1.0, .yz = 0.0,
                                               .zz = 1.0};
  lua_push_sym_tensor2(L, &I_sym);
  lua_setfield(L, -2, "unit");
  lua_pop(L, 1);

  return 0;
}

static int lua_register_mpi(lua_State* L)
{
  lua_newtable(L);
  lua_setglobal(L, "mpi");

  // Register the communicator class.
  lua_register_record_type(L, "mpi.comm", "An MPI communicator.", mpi_comm_funcs, mpi_comm_fields, mpi_comm_mm);

  // Now register MPI_COMM_WORLD and MPI_COMM_SELF objects.
  lua_getglobal(L, "mpi");
  lua_push_mpi_comm(L, MPI_COMM_WORLD);
  lua_setfield(L, -2, "COMM_WORLD");
  lua_push_mpi_comm(L, MPI_COMM_SELF);
  lua_setfield(L, -2, "COMM_SELF");

  return 0;
}

int lua_register_core_modules(lua_State* L)
{
  // Core types.
  lua_register_record_type(L, "complex", "A complex number.", complex_funcs, complex_fields, complex_mm);
  lua_register_record_type(L, "point", "A 3D point.", point_funcs, point_fields, point_mm);
  lua_register_record_type(L, "vector", "A 3D vector.", vector_funcs, vector_fields, vector_mm);
  lua_register_record_type(L, "tensor2", "A general rank 2 tensor.", tensor2_funcs, tensor2_fields, tensor2_mm);
  lua_register_record_type(L, "sym_tensor2", "A symmetric rank 2 tensor.", sym_tensor2_funcs, sym_tensor2_fields, sym_tensor2_mm);
  lua_register_array(L);
  lua_register_ndarray(L);
  lua_register_constants(L);
  lua_register_mpi(L);

  lua_register_class(L, "bbox", "A 3D bounding box.", bbox_funcs, bbox_methods);
  lua_register_class(L, "sp_func", "A function in 3D space.", sp_funcs, sp_methods);
  lua_register_class(L, "st_func", "A time-dependent function in 3D space.", st_funcs, st_methods);
  lua_register_class(L, "tagger", "An object that holds tags.", tagger_funcs, tagger_methods);
  lua_register_record_type(L, "point_cloud", "A point cloud in 3D space.", pc_funcs, pc_fields, pc_mm);

  // Register the options table.
  lua_register_options(L);

  // Utility functions.
  lua_register_util_funcs(L);

  return 0;
}

void lua_push_real(lua_State* L, real_t x)
{
  lua_pushnumber(L, (double)x);
}

bool lua_is_real(lua_State* L, int index)
{
  return lua_isnumber(L, index);
}

real_t lua_to_real(lua_State* L, int index)
{
  return (real_t)(lua_tonumber(L, index));
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

void lua_push_mpi_comm(lua_State* L, MPI_Comm comm)
{
  MPI_Comm* c = polymec_malloc(sizeof(MPI_Comm));
  *c = comm;
  lua_push_record(L, "mpi.comm", c, polymec_free);
}

bool lua_is_mpi_comm(lua_State* L, int index)
{
  return lua_is_record(L, index, "mpi.comm");
}

MPI_Comm lua_to_mpi_comm(lua_State* L, int index)
{
  return *((MPI_Comm*)lua_to_record(L, index, "mpi.comm"));
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

void lua_push_tagger(lua_State* L, tagger_t* t)
{
  lua_push_object(L, "tagger", t, NULL);
}

bool lua_is_tagger(lua_State* L, int index)
{
  return lua_is_object(L, index, "tagger");
}

tagger_t* lua_to_tagger(lua_State* L, int index)
{
  return (tagger_t*)lua_to_object(L, index, "tagger");
}

void lua_push_point_cloud(lua_State* L, point_cloud_t* c)
{
  lua_push_record(L, "point_cloud", c, DTOR(point_cloud_free));
}

bool lua_is_point_cloud(lua_State* L, int index)
{
  return lua_is_record(L, index, "point_cloud");
}

point_cloud_t* lua_to_point_cloud(lua_State* L, int index)
{
  return (point_cloud_t*)lua_to_record(L, index, "point_cloud");
}

