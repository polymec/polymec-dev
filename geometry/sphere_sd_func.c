// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/sphere_sd_func.h"

typedef struct
{
  point_t x;
  real_t r;
  normal_orient_t orient;
} sphere_t;

static real_t sphere_value(void* ctx, point_t* x)
{
  sphere_t* s = ctx;
  real_t sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  return sign * (point_distance(x, &s->x) - s->r);
}

static void sphere_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  sphere_t* s = ctx;
  real_t sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  real_t D = point_distance(x, &s->x);
  grad->x = sign * (x->x - s->x.x) / (D + 1e-14);
  grad->y = sign * (x->y - s->x.y) / (D + 1e-14);
  grad->z = sign * (x->z - s->x.z) / (D + 1e-14);
}

sd_func_t* sphere_sd_func_new(point_t* x, real_t r, normal_orient_t normal_orientation)
{
  // Set up a sphere signed distance function.
  sphere_t* s = polymec_malloc(sizeof(sphere_t));
  point_copy(&s->x, x);
  s->r = r;
  s->orient = normal_orientation;

  char sphere_str[1024];
  sprintf(sphere_str, "Sphere (x = (%g %g %g), r = %g)",
          x->x, x->y, x->z, r);
  sd_func_vtable vtable = {.value = sphere_value,
                           .eval_grad = sphere_eval_grad,
                           .dtor = polymec_free};
  return sd_func_new(sphere_str, s, vtable);
}

