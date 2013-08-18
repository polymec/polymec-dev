// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "geometry/sphere.h"

typedef struct
{
  point_t x;
  double r;
  normal_orient_t orient;
} sphere_t;

static void sphere_eval(void* ctx, point_t* x, double* result)
{
  sphere_t* s = ctx;
  double sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  result[0] = sign * (point_distance(x, &s->x) - s->r);
}

static void sphere_eval_gradient(void* ctx, point_t* x, double* result)
{
  sphere_t* s = ctx;
  double sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  double D = point_distance(x, &s->x);
  result[0] = sign * (x->x - s->x.x) / (D + 1e-14);
  result[1] = sign * (x->y - s->x.y) / (D + 1e-14);
  result[2] = sign * (x->z - s->x.z) / (D + 1e-14);
}

sp_func_t* sphere_new(point_t* x, double r, normal_orient_t normal_orientation)
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

