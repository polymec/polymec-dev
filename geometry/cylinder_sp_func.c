// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/cylinder_sp_func.h"

typedef struct
{
  vector_t d;
  point_t x;
  real_t r;
  normal_orient_t orient;
} cyl_t;

static void cyl_eval(void* ctx, point_t* x, real_t* result)
{
  cyl_t* c = ctx;
  real_t sign = (c->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  real_t r2 = (x->x - c->x.x)*(x->x - c->x.x) + 
              (x->y - c->x.y)*(x->y - c->x.y);
  result[0] = sign * (sqrt(r2) - c->r);
}

static void cyl_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  cyl_t* c = ctx;
  real_t r2 = (x->x - c->x.x)*(x->x - c->x.x) + 
              (x->y - c->x.y)*(x->y - c->x.y);
  real_t D = ABS(sqrt(r2) - c->r);
  if (reals_equal(D, 0.0))
    result[0] = result[1] = result[2] = 0.0;
  else
  {
    result[0] = (x->x - c->x.x) / D;
    result[1] = (x->y - c->x.y) / D;
    result[2] = 0.0;
  }
}

sp_func_t* cylinder_sp_func_new(point_t* x, real_t r, normal_orient_t normal_orientation)
{
  // Set up a cylinder signed distance function.
  cyl_t* c = polymec_malloc(sizeof(cyl_t));
  point_copy(&c->x, x);
  c->r = r;
  c->orient = normal_orientation;

  char cyl_str[1024];
  sprintf(cyl_str, "Cylinder (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_func_vtable vtable = {.eval = cyl_eval, .dtor = polymec_free};
  sp_func_t* cyl = sp_func_new(cyl_str, c, vtable, SP_FUNC_HETEROGENEOUS, 1);

  // Register the gradient function.
  char cyl_grad_str[1024];
  sprintf(cyl_grad_str, "Cylinder gradient (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_func_vtable vtable_g = {.eval = cyl_eval_gradient}; // Notice no dtor.
  sp_func_t* cyl_grad = sp_func_new(cyl_grad_str, c, vtable_g, SP_FUNC_HETEROGENEOUS, 3);
  sp_func_register_deriv(cyl, 1, cyl_grad);

  return cyl;
}

