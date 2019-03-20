// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/cylinder_sd_func.h"

typedef struct
{
  vector_t d;
  point_t x;
  real_t r;
  normal_orient_t orient;
} cyl_t;

static real_t cyl_value(void* ctx, point_t* x)
{
  cyl_t* c = ctx;
  real_t sign = (c->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  real_t r2 = (x->x - c->x.x)*(x->x - c->x.x) +
              (x->y - c->x.y)*(x->y - c->x.y);
  return sign * (sqrt(r2) - c->r);
}

static void cyl_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  cyl_t* c = ctx;
  real_t r2 = (x->x - c->x.x)*(x->x - c->x.x) +
              (x->y - c->x.y)*(x->y - c->x.y);
  real_t D = ABS(sqrt(r2) - c->r);
  if (reals_equal(D, 0.0))
    grad->x = grad->y = grad->z = 0.0;
  else
  {
    grad->x = (x->x - c->x.x) / D;
    grad->y = (x->y - c->x.y) / D;
    grad->z = 0.0;
  }
}

sd_func_t* cylinder_sd_func_new(point_t* x, real_t r, normal_orient_t normal_orientation)
{
  // Set up a cylinder signed distance function.
  cyl_t* c = polymec_malloc(sizeof(cyl_t));
  point_copy(&c->x, x);
  c->r = r;
  c->orient = normal_orientation;

  char cyl_str[1024];
  sprintf(cyl_str, "Cylinder (x = (%g %g %g), r = %g)",
          x->x, x->y, x->z, r);
  sd_func_vtable vtable = {.value = cyl_value,
                           .eval_grad = cyl_eval_grad,
                           .dtor = polymec_free};
  return sd_func_new(cyl_str, c, vtable);
}

