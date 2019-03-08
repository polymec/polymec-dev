// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/timer.h"
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

static lua_class_field complex_fields[] = {
  {"real", z_real, NULL},
  {"imag", z_imag, NULL},
  {NULL, NULL, NULL}
};

static int z_add(lua_State* L)
{
  if (!lua_is_complex(L, 1) || (!lua_is_complex(L, 2)))
    return luaL_error(L, "Arguments must both be complex numbers.");
  complex_t z1 = lua_to_complex(L, 1);
  complex_t z2 = lua_to_complex(L, 2);
  lua_push_complex(L, z1 + z2);
  return 1;
}

static int z_sub(lua_State* L)
{
  if (!lua_is_complex(L, 1) || (!lua_is_complex(L, 2)))
    return luaL_error(L, "Arguments must both be complex numbers.");
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

static lua_class_method complex_methods[] = {
  {"__add", z_add, NULL},
  {"__sub", z_sub, NULL},
  {"__mul", z_mul, NULL},
  {"__div", z_div, NULL},
  {"__unm", z_unm, NULL},
  {"__pow", z_pow, NULL},
  {"__tostring", z_tostring, NULL},
  {NULL, NULL, NULL}
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

static int p_random(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_is_rng(L, 1) || !lua_is_bbox(L, 2))
  {
    return luaL_error(L, "Arguments must be a rng and a bbox.");
  }

  rng_t* rng = lua_to_rng(L, 1);
  bbox_t* box = lua_to_bbox(L, 2);
  point_t* point = point_new(0.0, 0.0, 0.0);
  point_randomize(point, rng, box);
  lua_push_point(L, point);
  return 1;
}

static int p_distance(lua_State* L)
{
  point_t* p = lua_to_point(L, 1);
  if (p == NULL)
    return luaL_error(L, "Argument 1 must be a point.");
  point_t* q = lua_to_point(L, 2);
  if (q == NULL)
    return luaL_error(L, "Argument 2 must be a point.");
  lua_pushnumber(L, (double)point_distance(p, q));
  return 1;
}

static lua_module_function point_funcs[] = {
  {"new", p_new, "point.new(x, y, z) -> new 3D point."},
  {"random", p_random, "point.random(rng, bbox) -> new point with random coordinates within the given bounding box."},
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

static lua_class_field point_fields[] = {
  {"x", p_x, p_set_x},
  {"y", p_y, p_set_y},
  {"z", p_z, p_set_z},
  {NULL, NULL, NULL}
};

static int p_sub(lua_State* L)
{
  point_t* p1 = lua_to_point(L, 1);
  point_t* p2 = lua_to_point(L, 2);
  if ((p1 == NULL) || (p2 == NULL))
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

static lua_class_method point_methods[] = {
  {"__sub", p_sub, NULL},
  {"__len", p_len, NULL},
  {"__tostring", p_tostring, NULL},
  {NULL, NULL, NULL}
};

static int p2_new(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 2) ||
      !lua_isnumber(L, 1) || !lua_isnumber(L, 2))
  {
    return luaL_error(L, "Arguments must be x, y coordinates.");
  }

  real_t x = lua_to_real(L, 1);
  real_t y = lua_to_real(L, 2);
  point2_t* point = point2_new(x, y);
  lua_push_point2(L, point);
  return 1;
}

static int p2_distance(lua_State* L)
{
  point2_t* p = lua_to_point2(L, 1);
  point2_t* q = lua_to_point2(L, 2);
  if (q == NULL)
    return luaL_error(L, "Argument must be a point2.");
  lua_pushnumber(L, (double)point2_distance(p, q));
  return 1;
}

static lua_module_function point2_funcs[] = {
  {"new", p2_new, "point2.new(x, y, z) -> new 3D point."},
  {"distance", p2_distance, "point2.distance(x, y) -> |x - y|."},
  {NULL, NULL, NULL}
};

static int p2_x(lua_State* L)
{
  point2_t* p = lua_to_point2(L, 1);
  lua_pushnumber(L, (double)p->x);
  return 1;
}

static int p2_set_x(lua_State* L)
{
  point2_t* p = lua_to_point2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->x = lua_to_real(L, 2);
  return 0;
}

static int p2_y(lua_State* L)
{
  point2_t* p = lua_to_point2(L, 1);
  lua_pushnumber(L, (double)p->y);
  return 1;
}

static int p2_set_y(lua_State* L)
{
  point2_t* p = lua_to_point2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Point coordinates must be numbers.");
  p->y = lua_to_real(L, 2);
  return 0;
}

static lua_class_field point2_fields[] = {
  {"x", p2_x, p2_set_x},
  {"y", p2_y, p2_set_y},
  {NULL, NULL, NULL}
};

static int p2_len(lua_State* L)
{
  lua_pushinteger(L, 2);
  return 1;
}

static int p2_tostring(lua_State* L)
{
  point2_t* p = lua_to_point2(L, 1);
  lua_pushfstring(L, "point2 (%f, %f)", p->x, p->y);
  return 1;
}

static lua_class_method point2_methods[] = {
  {"__len", p2_len, NULL},
  {"__tostring", p2_tostring, NULL},
  {NULL, NULL, NULL}
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

static int v_random(lua_State* L)
{
  // Check the arguments.
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_is_rng(L, 1) || !lua_isnumber(L, 2))
    return luaL_error(L, "Arguments must be a rng and a (positive) magnitude.");

  rng_t* rng = lua_to_rng(L, 1);
  real_t mag = lua_to_real(L, 2);
  if (mag <= 0.0)
    return luaL_error(L, "Argument 2 must be a positive number.");
  vector_t* v = vector_new(0.0, 0.0, 0.0);
  vector_randomize(v, rng, mag);
  lua_push_vector(L, v);
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

static int v_cross(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (w == NULL)
    return luaL_error(L, "Argument must be a vector.");
  vector_t* vxw = vector_new(v->y*w->z - v->z*w->y,
                             v->z*w->x - v->x*w->z,
                             v->x*w->y - v->y*w->x);
  lua_push_vector(L, vxw);
  return 1;
}

static int v_dyad(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if (w == NULL)
    return luaL_error(L, "Argument must be a vector.");
  tensor2_t* vw = tensor2_new(v->x*w->x, v->x*w->y, v->x*w->z,
                              v->y*w->x, v->y*w->y, v->y*w->z,
                              v->z*w->x, v->z*w->y, v->z*w->z);
  lua_push_tensor2(L, vw);
  return 1;
}

static lua_module_function vector_funcs[] = {
  {"new", v_new, "vector.new(vx, vy, vz) -> new 3D vector."},
  {"random", v_random, "vector.random(rng, mag) -> new 3D vector with random components having the given magnitude."},
  {"dot", v_dot, "vector.dot(v1, v2) -> dot product of v1 with v2."},
  {"cross", v_cross, "vector.cross(v1, v2) -> cross product of v1 with v2."},
  {"dyad", v_dyad, "vector.dyad(v1, v2) -> dyadic product of v1 with v2."},
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

static int v_mag(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)vector_mag(v));
  return 1;
}

static int v_mag2(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  lua_pushnumber(L, (double)vector_dot(v, v));
  return 1;
}

static lua_class_field vector_fields[] = {
  {"x", v_x, v_set_x},
  {"y", v_y, v_set_y},
  {"z", v_z, v_set_z},
  {"mag", v_mag, NULL},
  {"mag2", v_mag2, NULL},
  {NULL, NULL, NULL}
};

static int v_add(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if ((v == NULL) || (w == NULL))
    luaL_error(L, "Arguments must both be vectors.");
  vector_t* sum = vector_new(v->x + w->x, v->y + w->y, v->z + w->z);
  lua_push_vector(L, sum);
  return 1;
}

static int v_sub(lua_State* L)
{
  vector_t* v = lua_to_vector(L, 1);
  vector_t* w = lua_to_vector(L, 2);
  if ((v == NULL) || (w == NULL))
    luaL_error(L, "Arguments must both be vectors.");
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
  lua_pushfstring(L, "vector(%f, %f, %f)", v->x, v->y, v->z);
  return 1;
}

static lua_class_method vector_methods[] = {
  {"__add", v_add, NULL},
  {"__sub", v_sub, NULL},
  {"__mul", v_mul, NULL},
  {"__div", v_div, NULL},
  {"__unm", v_unm, NULL},
  {"__len", v_len, NULL},
  {"__tostring", v_tostring, NULL},
  {"dot", v_dot, "v:dot(w) -> dot product of vectors v and w."},
  {"cross", v_cross, "v:cross(w) -> cross product of vectors v and w."},
  {"dyad", v_dyad, "v:dyad(w) -> dyadic product of vectors v and w."},
  {NULL, NULL, NULL}
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

static lua_class_field tensor2_fields[] = {
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
  if ((t1 == NULL) || (t2 == NULL))
    luaL_error(L, "Arguments must both be tensor2s.");
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
  if ((t1 == NULL) || (t2 == NULL))
    luaL_error(L, "Arguments must both be tensor2s.");
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
  lua_pushfstring(L, "tensor2((%f, %f, %f), (%f, %f, %f), (%f, %f, %f))",
                  t->xx, t->xy, t->xz, t->yx, t->yy, t->yz, t->zx, t->zy, t->zz);
  return 1;
}

static lua_class_method tensor2_methods[] = {
  {"__add", t2_add, NULL},
  {"__sub", t2_sub, NULL},
  {"__mul", t2_mul, NULL},
  {"__div", t2_div, NULL},
  {"__unm", t2_unm, NULL},
  {"__len", t2_len, NULL},
  {"__tostring", t2_tostring, NULL},
  {NULL, NULL, NULL}
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
  symtensor2_t* t = symtensor2_new(xx, xy, xz,
                                       yy, yz,
                                           zz);
  lua_push_symtensor2(L, t);
  return 1;
}

static lua_module_function symtensor2_funcs[] = {
  {"new", st2_new, "symtensor2.new(txx, txy, txz, tyy, tyz, tzz) -> new symmetric rank 2 tensor."},
  {NULL, NULL, NULL}
};

static int st2_xx(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->xx);
  return 1;
}

static int st2_set_xx(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xx = lua_to_real(L, 2);
  return 0;
}

static int st2_xy(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->xy);
  return 1;
}

static int st2_set_xy(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xy = lua_to_real(L, 2);
  return 0;
}

static int st2_xz(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->xz);
  return 1;
}

static int st2_set_xz(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xz = lua_to_real(L, 2);
  return 0;
}

static int st2_yx(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->xy);
  return 1;
}

static int st2_set_yx(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xy = lua_to_real(L, 2);
  return 0;
}

static int st2_yy(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->yy);
  return 1;
}

static int st2_set_yy(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yy = lua_to_real(L, 2);
  return 0;
}

static int st2_yz(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->yz);
  return 1;
}

static int st2_set_yz(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yz = lua_to_real(L, 2);
  return 0;
}

static int st2_zx(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->xz);
  return 1;
}

static int st2_set_zx(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->xz = lua_to_real(L, 2);
  return 0;
}

static int st2_zy(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->yz);
  return 1;
}

static int st2_set_zy(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->yz = lua_to_real(L, 2);
  return 0;
}

static int st2_zz(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushnumber(L, (double)t->zz);
  return 1;
}

static int st2_set_zz(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Tensor components must be numbers.");
  t->zz = lua_to_real(L, 2);
  return 0;
}

static lua_class_field symtensor2_fields[] = {
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
  symtensor2_t* t1 = lua_to_symtensor2(L, 1);
  symtensor2_t* t2 = lua_to_symtensor2(L, 2);
  if ((t1 == NULL) || (t2 == NULL))
    luaL_error(L, "Arguments must both be symtensor2s.");
  symtensor2_t* sum = symtensor2_new(t1->xx + t2->xx, t1->xy + t2->xy, t1->xz + t2->xz,
                                                      t1->yy + t2->yy, t1->yz + t2->yz,
                                                                       t1->zz + t2->zz);
  lua_push_symtensor2(L, sum);
  return 1;
}

static int st2_sub(lua_State* L)
{
  symtensor2_t* t1 = lua_to_symtensor2(L, 1);
  symtensor2_t* t2 = lua_to_symtensor2(L, 2);
  if ((t1 == NULL) || (t2 == NULL))
    luaL_error(L, "Arguments must both be symtensor2s.");
  symtensor2_t* diff = symtensor2_new(t1->xx - t2->xx, t1->xy - t2->xy, t1->xz - t2->xz,
                                                       t1->yy - t2->yy, t1->yz - t2->yz,
                                                                        t1->zz - t2->zz);
  lua_push_symtensor2(L, diff);
  return 1;
}

static int st2_mul(lua_State* L)
{
  if ((!lua_isnumber(L, 1) || !lua_is_symtensor2(L, 2)) &&
      (!lua_is_symtensor2(L, 1) || !lua_isnumber(L, 2)))
    luaL_error(L, "Arguments must be a symtensor2 and a number.");
  symtensor2_t* t = lua_to_symtensor2(L, (lua_isnumber(L, 1)) ? 2 : 1);
  real_t c = lua_to_real(L, (lua_isnumber(L, 1)) ? 1 : 2);
  symtensor2_t* t1 = symtensor2_new(c * t->xx, c * t->xy, c * t->xz,
                                               c * t->yy, c * t->yz,
                                                          c * t->zz);
  lua_push_symtensor2(L, t1);
  return 1;
}

static int st2_div(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  if (t == NULL)
    luaL_error(L, "Argument 1 must be a symtensor2.");
  if (!lua_isnumber(L, 2))
    luaL_error(L, "Argument 2 must be a number.");
  real_t c = lua_to_real(L, 2);
  symtensor2_t* t1 = symtensor2_new(t->xx/c, t->xy/c, t->xz/c,
                                             t->yy/c, t->yz/c,
                                                      t->zz/c);
  lua_push_symtensor2(L, t1);
  return 1;
}

static int st2_unm(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  symtensor2_t* t1 = symtensor2_new(-1.0 * t->xx, -1.0 * t->xy, -1.0 * t->xz,
                                                  -1.0 * t->yy, -1.0 * t->yz,
                                                                -1.0 * t->zz);
  lua_push_symtensor2(L, t1);
  return 1;
}

static int st2_len(lua_State* L)
{
  lua_pushinteger(L, 6);
  return 1;
}

static int st2_tostring(lua_State* L)
{
  symtensor2_t* t = lua_to_symtensor2(L, 1);
  lua_pushfstring(L, "symtensor2((%f, %f, %f), (%f, %f, %f), (%f, %f, %f))",
                  t->xx, t->xy, t->xz, t->xy, t->yy, t->yz, t->xz, t->yz, t->zz);
  return 1;
}

static lua_class_method symtensor2_methods[] = {
  {"__add", st2_add, NULL},
  {"__sub", st2_sub, NULL},
  {"__mul", st2_mul, NULL},
  {"__div", st2_div, NULL},
  {"__unm", st2_unm, NULL},
  {"__len", st2_len, NULL},
  {"__tostring", st2_tostring, NULL},
  {NULL, NULL, NULL}
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

static lua_class_field mpi_comm_fields[] = {
  {"size", mpi_comm_size, NULL},
  {"rank", mpi_comm_rank, NULL},
  {NULL, NULL, NULL}
};

static int mpi_comm_tostring(lua_State* L)
{
  if (lua_gettop(L) == 0)
    return luaL_error(L, "Method must be invoked with an mpi.comm.");
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

static lua_class_method mpi_comm_methods[] = {
  {"__tostring", mpi_comm_tostring, NULL},
  {NULL, NULL, NULL}
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

static lua_module_function bbox_funcs[] = {
  {"new", bb_new, "bbox.new{x1 = X1, x2 = X2, y1 = Y1, y2 = Y2, z1 = Z1, z2 = Z2} -> new bounding box."},
  {NULL, NULL, NULL}
};

static int bb_x1(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)bb->x1);
  return 1;
}

static int bb_set_x1(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding box limits must be numbers.");
  bb->x1 = lua_to_real(L, 2);
  return 0;
}

static int bb_x2(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)bb->x2);
  return 1;
}

static int bb_set_x2(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding box limits must be numbers.");
  bb->x2 = lua_to_real(L, 2);
  return 0;
}

static int bb_y1(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)bb->y1);
  return 1;
}

static int bb_set_y1(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding box limits must be numbers.");
  bb->y1 = lua_to_real(L, 2);
  return 0;
}

static int bb_y2(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)bb->y2);
  return 1;
}

static int bb_set_y2(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding box limits must be numbers.");
  bb->y2 = lua_to_real(L, 2);
  return 0;
}

static int bb_z1(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)bb->z1);
  return 1;
}

static int bb_set_z1(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding box limits must be numbers.");
  bb->z1 = lua_to_real(L, 2);
  return 0;
}

static int bb_z2(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  lua_pushnumber(L, (double)bb->z2);
  return 1;
}

static int bb_set_z2(lua_State* L)
{
  bbox_t* bb = lua_to_bbox(L, 1);
  if (!lua_isnumber(L, 2))
    return luaL_error(L, "Bounding box limits must be numbers.");
  bb->z2 = lua_to_real(L, 2);
  return 0;
}

static lua_class_field bbox_fields[] = {
  {"x1", bb_x1, bb_set_x1},
  {"x2", bb_x2, bb_set_x2},
  {"y1", bb_y1, bb_set_y1},
  {"y2", bb_y2, bb_set_y2},
  {"z1", bb_z1, bb_set_z1},
  {"z2", bb_z2, bb_set_z2},
  {NULL, NULL, NULL}
};

static int bb_contains(lua_State* L)
{
  bbox_t* b = lua_to_bbox(L, 1);
  if (b == NULL)
    return luaL_error(L, "Method must be invoked with a bbox.");
  point_t* p = lua_to_point(L, 2);
  if (p == NULL)
    return luaL_error(L, "Argument must be a point.");
  lua_pushboolean(L, bbox_contains(b, p));
  return 1;
}

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

static int sp_new(lua_State* L)
{
  // Check the argument.
  if (!lua_isfunction(L, 1))
    return luaL_error(L, "Argument 1 must be a Lua function.");
  if (!lua_isinteger(L, 2))
    return luaL_error(L, "Argument 2 must be a number of components.");
  lua_Integer n = lua_tointeger(L, 2);
  if (n < 0)
    return luaL_error(L, "Argument 2 must be a positive number.");
  sp_func_t* f = lua_as_sp_func(L, 1, (int)n);
  if (f != NULL)
  {
    lua_push_sp_func(L, f);
    return 1;
  }
  else
    return luaL_error(L, "Could not create an sp_func from arguments.");
}

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
  {"new", sp_new, "sp_func.new(f, n) -> returns an n-component spatial function implemented by the lua function f(x)."},
  {"constant", sp_constant, "sp_func.constant(F0) -> returns a function with the constant value F0."},
  {NULL, NULL, NULL}
};

static int sp_get_name(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  lua_pushstring(L, sp_func_name(f));
  return 1;
}

static int sp_set_name(lua_State* L)
{
  sp_func_t* f = lua_to_sp_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  sp_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static lua_class_field sp_fields[] = {
  {"name", sp_get_name, sp_set_name},
  {NULL, NULL, NULL}
};

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
  {"__gc", sp_gc, NULL},
  {"__len", sp_len, NULL},
  {"__call", sp_call, NULL},
  {"__tostring", sp_tostring, NULL},
  {NULL, NULL, NULL}
};

static int st_new(lua_State* L)
{
  // Check the argument.
  if (!lua_isfunction(L, 1))
    return luaL_error(L, "Argument 1 must be a Lua function.");
  if (!lua_isinteger(L, 2))
    return luaL_error(L, "Argument 2 must be a number of components.");
  lua_Integer n = lua_tointeger(L, 2);
  if (n < 0)
    return luaL_error(L, "Argument 2 must be a positive number.");
  st_func_t* f= lua_as_st_func(L, 1, (int)n);
  if (f != NULL)
  {
    lua_push_st_func(L, f);
    return 1;
  }
  else
    return luaL_error(L, "Could not create an st_func from arguments.");
}

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
  {"new", st_new, "st_func.new(f, n) -> returns an n-component space-time function implemented by the lua function f(x, t)."},
  {"constant", st_constant, "st_func.constant(F0) -> returns a function with the constant value F0."},
  {"from_sp_func", st_from_sp_func, "st_func.from_sp_func(f) -> returns a space-time function g(x, t) identical in form to the spatial function f(x)."},
  {NULL, NULL, NULL}
};

static int st_get_name(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  lua_pushstring(L, st_func_name(f));
  return 1;
}

static int st_set_name(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  if (!lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a string.");
  st_func_rename(f, lua_tostring(L, 2));
  return 0;
}

static lua_class_field st_fields[] = {
  {"name", st_get_name, st_set_name},
  {NULL, NULL, NULL}
};


static int st_freeze(lua_State* L)
{
  st_func_t* f = lua_to_st_func(L, 1);
  if (f == NULL)
    return luaL_error(L, "Method must be invoked with an st_func.");
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
  {"freeze", st_freeze, "f:freeze(t) -> freezes value of f at time t."},
  {"__gc", st_gc, NULL},
  {"__len", st_len, NULL},
  {"__call", st_call, NULL},
  {"__tostring", st_tostring, NULL},
  {NULL, NULL, NULL}
};

#ifdef _BSD_SOURCE
static int posix_rn(lua_State* L)
{
  if (!lua_isinteger(L, 1)
    luaL_error(L, "Argument must be a number of bytes.");
  lua_Integer size = lua_tointeger(L, 1);
  if (size <= 0)
    luaL_error(L, "Argument must be positive.");
  lua_push_rng(L, posix_rng_new((size_t)size));
  return 1;
}
#endif

static int rand_rn(lua_State* L)
{
  lua_push_rng(L, rand_rng_new());
  return 1;
}

static int host_rn(lua_State* L)
{
  lua_push_rng(L, host_rng_new());
  return 1;
}

static lua_module_function rng_funcs[] = {
#ifdef _BSD_SOURCE
  {"posix", posix_rn, "rng.posix(size) -> A POSIX (random()-based) random number generator with the given state size (in bytes)."},
#endif
  {"rand", rand_rn, "A rand()-based random number generator."},
  {"host", host_rn, "The best basic random number generator available."},
  {NULL, NULL, NULL}
};

static int rn_get_name(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  lua_pushstring(L, rng_name(r));
  return 1;
}

static int rn_get_min(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  lua_pushinteger(L, (lua_Integer)(rng_min(r)));
  return 1;
}

static int rn_get_max(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  lua_pushinteger(L, (lua_Integer)(rng_max(r)));
  return 1;
}

static lua_class_field rng_fields[] = {
  {"name", rn_get_name, NULL},
  {"min", rn_get_min, NULL},
  {"max", rn_get_max, NULL},
  {NULL, NULL, NULL}
};

static int rn_set_seed(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  if (r == NULL)
    luaL_error(L, "Method must be invoked with an rng.");
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Argument must be an unsigned 32-bit integer.");
  lua_Integer seed = lua_tointeger(L, 2);
  if (seed < 0)
    luaL_error(L, "Seed must be positive.");
  rng_set_seed(r, (uint32_t)seed);
  return 0;
}

static int rn_get(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  if (r == NULL)
    luaL_error(L, "Method must be invoked with an rng.");
  lua_pushinteger(L, (lua_Integer)(rng_get(r)));
  return 1;
}

static int rn_uniform(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  if (r == NULL)
    luaL_error(L, "Method must be invoked with an rng.");
  lua_push_real(L, rng_uniform(r));
  return 1;
}

static int rn_uniform_positive(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  if (r == NULL)
    luaL_error(L, "Method must be invoked with an rng.");
  lua_push_real(L, rng_uniform_positive(r));
  return 1;
}

static int rn_uniform_int(lua_State* L)
{
  rng_t* r = lua_to_rng(L, 1);
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Argument must be an unsigned 32-bit integer.");
  lua_Integer n = lua_tointeger(L, 2);
  if (n < 0)
    luaL_error(L, "n must be positive.");
  lua_pushinteger(L, (int)(rng_uniform_int(r, (uint32_t)n)));
  return 1;
}

static int rn_tostring(lua_State* L)
{
  rng_t* rng = lua_to_rng(L, 1);
  lua_pushfstring(L, "rng (%s)", rng_name(rng));
  return 1;
}

static lua_class_method rng_methods[] = {
  {"set_seed", rn_set_seed, "rng:set_seed(seed) -> Sets the seed for the random number generator."},
  {"get", rn_get, "rng:get() -> Generates and returns a random integer."},
  {"uniform", rn_uniform, "rng:uniform() -> Generates and returns a random floating point number in [0,1]."},
  {"uniform_positive", rn_uniform_positive, "rng:uniform_positive() -> Generates and returns a positive random floating point number in (0,1)."},
  {"uniform_int", rn_uniform_int, "rng:uniform_int(n) -> Generates and returns a random integer in the range [0, n-1]."},
  {"__tostring", rn_tostring, NULL},
  {NULL, NULL, NULL}
};

// FIXME: No adj_graph constructors yet!
static lua_module_function adj_funcs[] = {
  {NULL, NULL, NULL}
};

static int adj_get_comm(lua_State* L)
{
  adj_graph_t* g = lua_to_adj_graph(L, 1);
  lua_push_mpi_comm(L, adj_graph_comm(g));
  return 1;
}

static int adj_get_num_vertices(lua_State* L)
{
  adj_graph_t* g = lua_to_adj_graph(L, 1);
  lua_pushinteger(L, (lua_Integer)(adj_graph_num_vertices(g)));
  return 1;
}

static lua_class_field adj_fields[] = {
  {"comm", adj_get_comm, NULL},
  {"num_vertices", adj_get_num_vertices, NULL},
  {NULL, NULL, NULL}
};

static int adj_num_edges(lua_State* L)
{
  adj_graph_t* g = lua_to_adj_graph(L, 1);
  if (g == NULL)
    luaL_error(L, "Method must be invoked with an adj_graph.");
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Argument must be a vertex index.");
  int vertex = (int)lua_tointeger(L, 2);
  if (vertex < 0)
    luaL_error(L, "Vertex must be non-negative.");
  if (vertex >= adj_graph_num_vertices(g))
    luaL_error(L, "Invalid vertex: %d", vertex);
  lua_pushinteger(L, (lua_Integer)adj_graph_num_edges(g, vertex));
  return 0;
}

static int adj_edges(lua_State* L)
{
  adj_graph_t* g = lua_to_adj_graph(L, 1);
  if (g == NULL)
    luaL_error(L, "Method must be invoked with an adj_graph.");
  if (!lua_isinteger(L, 2))
    luaL_error(L, "Argument must be a vertex index.");
  int vertex = (int)lua_tointeger(L, 2);
  if (vertex < 0)
    luaL_error(L, "Vertex must be non-negative.");
  if (vertex >= adj_graph_num_vertices(g))
    luaL_error(L, "Invalid vertex: %d", vertex);
  size_t num_edges = adj_graph_num_edges(g, vertex);
  int* edges = adj_graph_edges(g, vertex);
  lua_newtable(L);
  for (size_t i = 0; i < num_edges; ++i)
  {
    lua_pushinteger(L, (lua_Integer)(edges[i]));
    lua_seti(L, -2, (int)i);
  }
  return 1;
}

static int adj_tostring(lua_State* L)
{
  adj_graph_t* g = lua_to_adj_graph(L, 1);
  lua_pushfstring(L, "adj_graph (%d vertices)", adj_graph_num_vertices(g));
  return 1;
}

static lua_class_method adj_methods[] = {
  {"num_edges", adj_num_edges, "graph:num_edges(vertex) -> Returns the number of edges for the given vertex."},
  {"edges", adj_edges, "graph:edges(vertex) -> Returns the set of edges for the given vertex."},
  {"__tostring", adj_tostring, NULL},
  {NULL, NULL, NULL}
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

// Here's a reals_equal(a, b) function just like our C version.
static int lua_reals_equal(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_is_real(L, 1) || !lua_is_real(L, 2))
    luaL_error(L, "reals_equal() accepts exactly two real arguments.");
  real_t a = lua_to_real(L, 1);
  real_t b = lua_to_real(L, 2);
  lua_pushboolean(L, reals_equal(a, b));
  return 1;
}

// And here's reals_nearly_equal(a, b, epsilon).
static int lua_reals_nearly_equal(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 3) || !lua_is_real(L, 1) || !lua_is_real(L, 2) || !lua_is_real(L, 3))
    luaL_error(L, "reals_equal() accepts exactly three real arguments.");
  real_t a = lua_to_real(L, 1);
  real_t b = lua_to_real(L, 2);
  real_t eps = lua_to_real(L, 3);
  lua_pushboolean(L, reals_nearly_equal(a, b, eps));
  return 1;
}

// This lets you set the tolerance for reals_equal() from Lua.
static int lua_set_real_epsilon(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) || !lua_is_real(L, 1))
    luaL_error(L, "set_real_epsilon() accepts exactly one real arguments.");
  real_t eps = lua_to_real(L, 1);
  set_real_epsilon(eps);
  return 0;
}

// Here's a dir() function just like Python's.
static int lua_dir(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 0) && (num_args != 1))
    luaL_error(L, "dir() accepts exactly zero or one arguments.");

  if (num_args == 0)
    lua_getglobal(L, "_G");

  int table_index = 0;
  if (lua_istable(L, 1))
    table_index = 1;
  else if (lua_isuserdata(L, 1) && (lua_getmetatable(L, 1) != 0))
    table_index = 2;
  else if (lua_isnil(L, 1))
  {
    lua_pushnil(L);
    return 1;
  }
  else
  {
    // Return an empty table.
    lua_newtable(L);
    return 1;
  }

  // Go through the keys of the table and push them into a new table.
  lua_newtable(L);
  int new_table_index = lua_gettop(L);
  int index = 0;
  lua_pushnil(L);
  while (lua_next(L, table_index))
  {
    ++index;

    // Key is at index -2, value is at -1.
    // First, replace the value with a copy of the key.
    lua_pop(L, 1);
    lua_pushvalue(L, -1);

    // Now add this key to our new table. This pops the key copy
    // off the stack.
    lua_rawseti(L, new_table_index, index);
  }

  // If this is a table *AND* it has a metatable, add the keys from the
  // metatable, too.
  if ((table_index == 1) && (lua_getmetatable(L, 1) != 0))
  {
    lua_pushnil(L);
    while (lua_next(L, -2))
    {
      ++index;

      // Key is at index -2, value is at -1.
      // First, replace the value with a copy of the key.
      lua_pop(L, 1);
      lua_pushvalue(L, -1);

      // Now add this key to our new table. This pops the key copy
      // off the stack.
      lua_rawseti(L, new_table_index, index);
    }
    lua_pop(L, 1); // pop the metatable
  }

  // Sort the table.
  lua_getglobal(L, "table");
  lua_getfield(L, -1, "sort");
  lua_pushvalue(L, new_table_index);
  lua_call(L, 1, 0);

  // Clean up the stack. The "table" module is at the top of the stack,
  // and our new table is underneath it, so we rotate once in the direction
  // of the top of the stack to bring our result to the top. This essentially
  // swaps the two indices.
  lua_rotate(L, lua_gettop(L)-1, 1);

  // Now just get rid of all the garbage in between our result and the
  // original arguments.
  while (lua_gettop(L) > (num_args + 1))
    lua_remove(L, -2);

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

  // Retrieve the object for the argument.
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
  // reals_equal(a, b) returns true if a == b to within polymec's tolerance,
  // false if not.
  lua_pushcfunction(L, lua_reals_equal);
  lua_setglobal(L, "reals_equal");
  lua_getglobal(L, "reals_equal");
  lua_set_docstring(L, -1, "reals_equal(a, b) -> Returns true if a == b for a, b real, false if not.");
  lua_pop(L, 1);

  // reals_nearly_equal(a, b, eps) returns true if a == b to within eps,
  // false if not.
  lua_pushcfunction(L, lua_reals_nearly_equal);
  lua_setglobal(L, "reals_nearly_equal");
  lua_getglobal(L, "reals_nearly_equal");
  lua_set_docstring(L, -1, "reals_nearly_equal(a, b, tolerance) -> Returns true if |a-b| < tolerance for a, b, tolerance real, false if not.");
  lua_pop(L, 1);

  lua_pushcfunction(L, lua_set_real_epsilon);
  lua_setglobal(L, "set_real_epsilon");
  lua_getglobal(L, "set_real_epsilon");
  lua_set_docstring(L, -1, "set_real_epsilon(epsilon) -> Sets the floating point tolerance.");
  lua_pop(L, 1);

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
  tensor2_t* I = tensor2_new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  lua_push_tensor2(L, I);
  lua_setfield(L, -2, "unit");
  lua_pop(L, 1);

  lua_getglobal(L, "symtensor2");
  symtensor2_t I_sym = {.xx = 1.0, .xy = 0.0, .xz = 0.0,
                                    .yy = 1.0, .yz = 0.0,
                                               .zz = 1.0};
  lua_push_symtensor2(L, &I_sym);
  lua_setfield(L, -2, "unit");
  lua_pop(L, 1);

  return 0;
}

static int mpi_get_comm_world(lua_State* L)
{
  lua_push_mpi_comm(L, MPI_COMM_WORLD);
  return 1;
}

static int mpi_get_comm_self(lua_State* L)
{
  lua_push_mpi_comm(L, MPI_COMM_SELF);
  return 1;
}

static lua_module_field mpi_fields[] =
{
  {"COMM_WORLD", mpi_get_comm_world, NULL},
  {"COMM_SELF", mpi_get_comm_self, NULL},
  {NULL, NULL, NULL}
};

static lua_module_function mpi_funcs[] = {
  {NULL, NULL, NULL}
};

static int lua_register_mpi(lua_State* L)
{
  lua_register_module(L, "mpi", "Message Passing Interface (MPI) tools.",
                      mpi_fields, mpi_funcs);

  // Register the communicator class.
  lua_register_class(L, "mpi.comm", "An MPI communicator.",
                     mpi_comm_funcs, mpi_comm_fields, mpi_comm_methods,
                     polymec_free);

  return 0;
}


static int l_get_level(lua_State* L)
{
  log_level_t level = log_level();
  if (level == LOG_NONE)
    lua_pushstring(L, "none");
  else if (level == LOG_URGENT)
    lua_pushstring(L, "urgent");
  else if (level == LOG_INFO)
    lua_pushstring(L, "info");
  else if (level == LOG_DETAIL)
    lua_pushstring(L, "detail");
  else // if (level == LOG_DEBUG)
    lua_pushstring(L, "debug");
  return 1;
}

static int l_set_level(lua_State* L)
{
  if (!lua_isstring(L, 2))
    return luaL_error(L, "log level must be 'none', 'urgent', 'info', 'detail', or 'debug'.");
  const char* level_str = lua_tostring(L, 2);
  log_level_t level;
  if (strcmp(level_str, "none") == 0)
    level = LOG_NONE;
  else if (strcmp(level_str, "urgent") == 0)
    level = LOG_URGENT;
  else if (strcmp(level_str, "info") == 0)
    level = LOG_INFO;
  else if (strcmp(level_str, "detail") == 0)
    level = LOG_DETAIL;
  else if (strcmp(level_str, "debug") == 0)
    level = LOG_DEBUG;
  else
    return luaL_error(L, "log level must be 'none', 'urgent', 'info', 'detail', or 'debug'.");
  set_log_level(level);
  return 0;
}

static int l_get_mode(lua_State* L)
{
  log_mode_t mode = log_mode();
  if (mode == LOG_TO_SINGLE_RANK)
    lua_pushstring(L, "single_rank");
  else // if (mode == LOG_TO_ALL_RANKS)
    lua_pushstring(L, "all_ranks");
  return 1;
}

static int l_set_mode(lua_State* L)
{
  if (!lua_isstring(L, 2))
    return luaL_error(L, "log mode must be 'single_rank' or 'all_ranks'.");
  const char* mode_str = lua_tostring(L, 2);
  log_mode_t mode;
  if (strcmp(mode_str, "single_rank") == 0)
    mode = LOG_TO_SINGLE_RANK;
  else if (strcmp(mode_str, "all_ranks") == 0)
    mode = LOG_TO_ALL_RANKS;
  else
    return luaL_error(L, "log mode must be 'single_rank' or 'all_ranks'.");
  set_log_mode(mode);
  return 0;
}

static int l_get_rank(lua_State* L)
{
  int rank = log_mpi_rank(log_level());
  if (rank != -1)
    lua_pushinteger(L, rank);
  else
    lua_pushnil(L);
  return 1;
}

static int l_set_rank(lua_State* L)
{
  log_mode_t mode = log_mode();
  if (mode == LOG_TO_ALL_RANKS)
    return luaL_error(L, "cannot specify a log rank in 'all_ranks' log mode.");
  else
  {
    if (!lua_isinteger(L, 2))
      return luaL_error(L, "rank must be an integer.");
    int rank = (int)lua_tointeger(L, 2);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if ((rank < 0) || (rank >= nprocs))
      return luaL_error(L, "rank must be an integer between 0 and %d.", nprocs);
    set_log_mpi_rank(log_level(), rank);
  }
  return 0;
}

static int l_get_stream(lua_State* L)
{
  FILE* stream = log_stream(log_level());
  if (stream != NULL)
  {
    luaL_Stream* s = lua_newuserdata(L, sizeof(luaL_Stream));
    s->f = stream;
    s->closef = NULL; // not responsible for closing!
    luaL_setmetatable(L, LUA_FILEHANDLE);
  }
  else
    lua_pushnil(L);
  return 1;
}

static int l_set_stream(lua_State* L)
{
  if (lua_isnil(L, 2))
    set_log_stream(log_level(), NULL);
  else
  {
    luaL_Stream* stream = luaL_checkudata(L, 2, LUA_FILEHANDLE);
    if (stream == NULL)
      return luaL_error(L, "stream must be a file handle.");
    set_log_stream(log_level(), stream->f);
    stream->closef = NULL; // we assume responsibility. FIXME: Cannot guarantee is fclose!
  }
  return 0;
}

static lua_module_field logging_fields[] =
{
  {"level", l_get_level, l_set_level},
  {"mode", l_get_mode, l_set_mode},
  {"rank", l_get_rank, l_set_rank},
  {"stream", l_get_stream, l_set_stream},
  {NULL, NULL, NULL}
};

static int lua_register_logging(lua_State* L)
{
  lua_register_module(L, "logging", "Settings for logging.",
                      logging_fields, NULL);
  return 0;
}

static int timers_get_enabled(lua_State* L)
{
  options_t* opts = options_argv();
  char* timers_on = options_value(opts, "timers");
  lua_pushboolean(L, (timers_on != NULL) && string_as_boolean(timers_on));
  return 1;
}

static int timers_get_file(lua_State* L)
{
  options_t* opts = options_argv();
  char* timers_on = options_value(opts, "timers");
  if ((timers_on != NULL) && string_as_boolean(timers_on))
    lua_pushstring(L, polymec_timer_file());
  else
    lua_pushnil(L);
  return 1;
}

static int timers_set_file(lua_State* L)
{
  if (lua_isnil(L, 2))
  {
    // Set to default value.
    polymec_set_timer_file("timer_report.txt");
  }
  else if (!lua_isstring(L, 2))
    return luaL_error(L, "timers.file must be the name of a file.");
  else
    polymec_set_timer_file(lua_tostring(L, 2));
  return 0;
}

static lua_module_field timers_fields[] =
{
  {"enabled", timers_get_enabled, NULL},
  {"file", timers_get_file, timers_set_file},
  {NULL, NULL, NULL}
};

static int lua_register_timers(lua_State* L)
{
  lua_register_module(L, "timers", "Settings for timers.",
                      timers_fields, NULL);
  return 0;
}

static int dl_get_paths(lua_State* L)
{
  int pos = 0;
  const char* path;
  lua_newtable(L);
  while (polymec_next_dl_path(&pos, &path))
  {
    lua_pushstring(L, path);
    lua_seti(L, -2, pos);
  }
  return 1;
}

static lua_module_field dl_fields[] =
{
  {"paths", dl_get_paths, NULL},
  {NULL, NULL, NULL}
};

static int dl_add_path(lua_State* L)
{
  if (!lua_isstring(L, 1))
    luaL_error(L, "Path must be an absolute path.");
  const char* path = lua_tostring(L, 1);
  if (path[0] != '/')
    luaL_error(L, "Path must be an absolute path.");
  polymec_add_dl_path(path);
  return 0;
}

static lua_module_function dl_funcs[] = {
  {"add_path", dl_add_path, "dl.add_path(path) -> Adds the given path to the list of paths searched for dynamically loadable libraries."},
  {NULL, NULL, NULL}
};

static int lua_register_dl(lua_State* L)
{
  lua_register_module(L, "dl", "Settings for dynamic load libraries.",
                      dl_fields, dl_funcs);
  return 0;
}

int lua_register_core_modules(lua_State* L)
{
  // Core types.
  lua_register_class(L, "complex", "A complex number.", complex_funcs, complex_fields, complex_methods, polymec_free);
  lua_register_class(L, "point", "A 3D point.", point_funcs, point_fields, point_methods, NULL);
  lua_register_class(L, "point2", "A 2D point.", point2_funcs, point2_fields, point2_methods, NULL);
  lua_register_class(L, "vector", "A 3D vector.", vector_funcs, vector_fields, vector_methods, NULL);
  lua_register_class(L, "tensor2", "A general rank 2 tensor.", tensor2_funcs, tensor2_fields, tensor2_methods, NULL);
  lua_register_class(L, "symtensor2", "A symmetric rank 2 tensor.", symtensor2_funcs, symtensor2_fields, symtensor2_methods, NULL);
  lua_register_array(L);
  lua_register_ndarray(L);
  lua_register_constants(L);
  lua_register_mpi(L);
  lua_register_logging(L);
  lua_register_timers(L);
  lua_register_dl(L);

  lua_register_class(L, "bbox", "A 3D bounding box.", bbox_funcs, bbox_fields, bbox_methods, NULL);
  lua_register_class(L, "sp_func", "A function in 3D space.", sp_funcs, sp_fields, sp_methods, NULL);
  lua_register_class(L, "st_func", "A time-dependent function in 3D space.", st_funcs, st_fields, st_methods, NULL);
  lua_register_class(L, "rng", "A random number generator.", rng_funcs, rng_fields, rng_methods, NULL);
  lua_register_class(L, "adj_graph", "An adjacency graph.", adj_funcs, adj_fields, adj_methods, NULL);

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
  lua_push_object(L, "complex", zz);
}

bool lua_is_complex(lua_State* L, int index)
{
  return lua_is_object(L, index, "complex");
}

complex_t lua_to_complex(lua_State* L, int index)
{
  return *((complex_t*)lua_to_object(L, index, "complex"));
}

void lua_push_point(lua_State* L, point_t* p)
{
  lua_push_object(L, "point", p);
}

bool lua_is_point(lua_State* L, int index)
{
  return lua_is_object(L, index, "point");
}

point_t* lua_to_point(lua_State* L, int index)
{
  return (point_t*)lua_to_object(L, index, "point");
}

void lua_push_point2(lua_State* L, point2_t* p)
{
  lua_push_object(L, "point2", p);
}

bool lua_is_point2(lua_State* L, int index)
{
  return lua_is_object(L, index, "point2");
}

point2_t* lua_to_point2(lua_State* L, int index)
{
  return (point2_t*)lua_to_object(L, index, "point2");
}

void lua_push_vector(lua_State* L, vector_t* v)
{
  lua_push_object(L, "vector", v);
}

bool lua_is_vector(lua_State* L, int index)
{
  return lua_is_object(L, index, "vector");
}

vector_t* lua_to_vector(lua_State* L, int index)
{
  return (vector_t*)lua_to_object(L, index, "vector");
}

void lua_push_tensor2(lua_State* L, tensor2_t* t)
{
  lua_push_object(L, "tensor2", t);
}

bool lua_is_tensor2(lua_State* L, int index)
{
  return lua_is_object(L, index, "tensor2");
}

tensor2_t* lua_to_tensor2(lua_State* L, int index)
{
  return (tensor2_t*)lua_to_object(L, index, "tensor2");
}

void lua_push_symtensor2(lua_State* L, symtensor2_t* t)
{
  lua_push_object(L, "symtensor2", t);
}

bool lua_is_symtensor2(lua_State* L, int index)
{
  return lua_is_object(L, index, "symtensor2");
}

symtensor2_t* lua_to_symtensor2(lua_State* L, int index)
{
  return (symtensor2_t*)lua_to_object(L, index, "symtensor2");
}

void lua_push_mpi_comm(lua_State* L, MPI_Comm comm)
{
  MPI_Comm* c = polymec_malloc(sizeof(MPI_Comm));
  *c = comm;
  lua_push_object(L, "mpi.comm", c);
}

bool lua_is_mpi_comm(lua_State* L, int index)
{
  return lua_is_object(L, index, "mpi.comm");
}

MPI_Comm lua_to_mpi_comm(lua_State* L, int index)
{
  return *((MPI_Comm*)lua_to_object(L, index, "mpi.comm"));
}

void lua_push_bbox(lua_State* L, bbox_t* b)
{
  lua_push_object(L, "bbox", b);
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
  lua_push_object(L, "sp_func", f);
}

bool lua_is_sp_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "sp_func");
}

sp_func_t* lua_to_sp_func(lua_State* L, int index)
{
  return (sp_func_t*)lua_to_object(L, index, "sp_func");
}

typedef struct
{
  lua_State* L;
  int num_comp;
} lua_sp_t;

static void lua_sp_eval(void* context, point_t* x, real_t* F)
{
  // Fetch our Lua table from the registry.
  lua_sp_t* sp = context;
  lua_pushlightuserdata(sp->L, sp);
  lua_gettable(sp->L, LUA_REGISTRYINDEX);

  // Call the function with x as an argument.
  lua_push_point(sp->L, x);
  lua_call(sp->L, 1, sp->num_comp);

  for (int c = 0; c < sp->num_comp; ++c)
  {
    if (!lua_isnumber(sp->L, -1))
      luaL_error(sp->L, "sp_func result is not a number.");
    F[c] = lua_to_real(sp->L, -1);
    lua_pop(sp->L, 1);
  }
}

static void lua_sp_dtor(void* context)
{
  lua_sp_t* sp = context;
  lua_pushnil(sp->L);
  lua_settable(sp->L, LUA_REGISTRYINDEX);
  polymec_free(sp);
}

sp_func_t* lua_as_sp_func(lua_State* L, int index, int num_comp)
{
  if (lua_isfunction(L, index))
  {
    lua_sp_t* sp = polymec_malloc(sizeof(lua_sp_t));
    sp->L = L;
    sp->num_comp = num_comp;
    sp_func_vtable vtable = {.eval = lua_sp_eval, .dtor = lua_sp_dtor};
    lua_pushlightuserdata(L, sp);
    lua_pushvalue(L, index);
    lua_rawset(L, LUA_REGISTRYINDEX);
    return sp_func_new("Lua sp_func", sp, vtable, SP_FUNC_HETEROGENEOUS, num_comp);
  }
  else
    return NULL;
}

void lua_push_st_func(lua_State* L, st_func_t* f)
{
  lua_push_object(L, "st_func", f);
}

bool lua_is_st_func(lua_State* L, int index)
{
  return lua_is_object(L, index, "st_func");
}

st_func_t* lua_to_st_func(lua_State* L, int index)
{
  return (st_func_t*)lua_to_object(L, index, "st_func");
}

static void lua_st_eval(void* context, point_t* x, real_t t, real_t* F)
{
  // Fetch our Lua table from the registry.
  lua_sp_t* st = context;
  lua_pushlightuserdata(st->L, st);
  lua_gettable(st->L, LUA_REGISTRYINDEX);

  // Call the function with x as an argument.
  lua_push_point(st->L, x);
  lua_push_real(st->L, t);
  lua_call(st->L, 2, st->num_comp);

  for (int c = 0; c < st->num_comp; ++c)
  {
    if (!lua_isnumber(st->L, -1))
      luaL_error(st->L, "st_func result is not a number.");
    F[c] = lua_to_real(st->L, -1);
    lua_pop(st->L, 1);
  }
}

st_func_t* lua_as_st_func(lua_State* L, int index, int num_comp)
{
  if (lua_isfunction(L, index))
  {
    lua_sp_t* st = polymec_malloc(sizeof(lua_sp_t));
    st->L = L;
    st->num_comp = num_comp;
    st_func_vtable vtable = {.eval = lua_st_eval, .dtor = lua_sp_dtor};
    lua_pushlightuserdata(L, st);
    lua_pushvalue(L, index);
    lua_rawset(L, LUA_REGISTRYINDEX);
    return st_func_new("Lua st_func", st, vtable,
                       ST_FUNC_HETEROGENEOUS, ST_FUNC_NONCONSTANT,
                       num_comp);
  }
  else
    return NULL;
}

void lua_push_rng(lua_State* L, rng_t* r)
{
  lua_push_object(L, "rng", r);
}

bool lua_is_rng(lua_State* L, int index)
{
  return lua_is_object(L, index, "rng");
}

rng_t* lua_to_rng(lua_State* L, int index)
{
  return (rng_t*)lua_to_object(L, index, "rng");
}

void lua_push_adj_graph(lua_State* L, adj_graph_t* g)
{
  lua_push_object(L, "adj_graph", g);
}

bool lua_is_adj_graph(lua_State* L, int index)
{
  return lua_is_object(L, index, "adj_graph");
}

adj_graph_t* lua_to_adj_graph(lua_State* L, int index)
{
  return (adj_graph_t*)lua_to_object(L, index, "adj_graph");
}

