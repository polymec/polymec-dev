// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/sphere.h"

typedef struct
{
  point_t x;
  real_t r;
  normal_orient_t orient;
} sphere_t;

static void sphere_eval(void* ctx, point_t* x, real_t* result)
{
  sphere_t* s = ctx;
  real_t sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  result[0] = sign * (point_distance(x, &s->x) - s->r);
}

static void sphere_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  sphere_t* s = ctx;
  real_t sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  real_t D = point_distance(x, &s->x);
  result[0] = sign * (x->x - s->x.x) / (D + 1e-14);
  result[1] = sign * (x->y - s->x.y) / (D + 1e-14);
  result[2] = sign * (x->z - s->x.z) / (D + 1e-14);
}

sp_func_t* sphere_new(point_t* x, real_t r, normal_orient_t normal_orientation)
{
  // Set up a sphere signed distance function.
  sphere_t* s = polymec_malloc(sizeof(sphere_t));
  point_copy(&s->x, x);
  s->r = r;
  s->orient = normal_orientation;

  char sphere_str[1024];
  sprintf(sphere_str, "Sphere (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_func_vtable vtable = {.eval = sphere_eval, .dtor = polymec_free};
  sp_func_t* sphere = sp_func_new(sphere_str, s, vtable, SP_FUNC_HETEROGENEOUS, 1);

  // Register the gradient function.
  char sphere_grad_str[1024];
  sprintf(sphere_grad_str, "Sphere gradient (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_func_vtable vtable_g = {.eval = sphere_eval_gradient}; // Notice no dtor.
  sp_func_t* sphere_grad = sp_func_new(sphere_grad_str, s, vtable_g, SP_FUNC_HETEROGENEOUS, 3);
  sp_func_register_deriv(sphere, 1, sphere_grad);

  return sphere;
}

