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
  sphere_t* s = malloc(sizeof(sphere_t));
  point_copy(&s->x, x);
  s->r = r;
  s->orient = normal_orientation;

  char sphere_str[1024];
  sprintf(sphere_str, "Sphere (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_vtable vtable = {.eval = sphere_eval, .dtor = free};
  sp_func_t* sphere = sp_func_new(sphere_str, s, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function.
  char sphere_grad_str[1024];
  sprintf(sphere_grad_str, "Sphere gradient (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_vtable vtable_g = {.eval = sphere_eval_gradient}; // Notice no dtor.
  sp_func_t* sphere_grad = sp_func_new(sphere_grad_str, s, vtable_g, SP_INHOMOGENEOUS, 3);
  sp_func_register_deriv(sphere, 1, sphere_grad);

  return sphere;
}

