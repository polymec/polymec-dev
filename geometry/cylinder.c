// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "geometry/cylinder.h"

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
  result[0] = sign * (rsqrt(r2) - c->r);
}

static void cyl_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  cyl_t* c = ctx;
  real_t r2 = (x->x - c->x.x)*(x->x - c->x.x) + 
              (x->y - c->x.y)*(x->y - c->x.y);
  real_t D = fabs(rsqrt(r2) - c->r);
  if (D == 0.0)
    result[0] = result[1] = result[2] = 0.0;
  else
  {
    result[0] = (x->x - c->x.x) / D;
    result[1] = (x->y - c->x.y) / D;
    result[2] = 0.0;
  }
}

sp_func_t* cylinder_new(point_t* x, real_t r, normal_orient_t normal_orientation)
{
  // Set up a cylinder signed distance function.
  cyl_t* c = malloc(sizeof(cyl_t));
  point_copy(&c->x, x);
  c->r = r;
  c->orient = normal_orientation;

  char cyl_str[1024];
  sprintf(cyl_str, "Cylinder (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_vtable vtable = {.eval = cyl_eval, .dtor = free};
  sp_func_t* cyl = sp_func_new(cyl_str, c, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function.
  char cyl_grad_str[1024];
  sprintf(cyl_grad_str, "Cylinder gradient (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_vtable vtable_g = {.eval = cyl_eval_gradient}; // Notice no dtor.
  sp_func_t* cyl_grad = sp_func_new(cyl_grad_str, c, vtable_g, SP_INHOMOGENEOUS, 3);
  sp_func_register_deriv(cyl, 1, cyl_grad);

  return cyl;
}

