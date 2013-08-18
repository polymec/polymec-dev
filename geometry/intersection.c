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

#include "geometry/intersection.h"

typedef struct
{
  sp_func_t** funcs;
  int num_funcs;
} inter_t;

// Destructor.
static void inter_free(void* ctx)
{
  inter_t* inter = ctx;
  free(inter->funcs);
  free(inter);
}

static void inter_eval(void* ctx, point_t* x, double* result)
{
  inter_t* inter = ctx;
  double maxval = -FLT_MAX;
  for (int i = 0; i < inter->num_funcs; ++i)
  {
    double ival;
    sp_func_eval(inter->funcs[i], x, &ival);
    maxval = MAX(maxval, ival);
  }
  result[0] = maxval;
}

static void inter_eval_gradient(void* ctx, point_t* x, double* result)
{
  inter_t* inter = ctx;
  double maxval = -FLT_MAX;
  int index = -1;
  for (int i = 0; i < inter->num_funcs; ++i)
  {
    double ival;
    sp_func_eval(inter->funcs[i], x, &ival);
    if (ival > maxval)
    {
      maxval = ival;
      index = i;
    }
  }
  ASSERT(sp_func_has_deriv(inter->funcs[index], 1));
  sp_func_eval_deriv(inter->funcs[index], 1, x, result);
}

sp_func_t* intersection_new(sp_func_t** surfaces, int num_surfaces)
{
  ASSERT(surfaces != NULL);
  ASSERT(num_surfaces > 1);

  inter_t* inter = malloc(sizeof(inter_t));
  inter->funcs = malloc(sizeof(sp_func_t*) * num_surfaces);
  inter->num_funcs = num_surfaces;
  bool has_grad = true;
  for (int i = 0; i < num_surfaces; ++i)
  {
    ASSERT(sp_func_num_comp(surfaces[i]) == 1);
    inter->funcs[i] = surfaces[i];
    if (!sp_func_has_deriv(surfaces[i], 1))
      has_grad = false;
  }

  char inter_str[1024];
  sprintf(inter_str, "Intersection"); // FIXME: Not very helpful.
  sp_vtable vtable = {.eval = inter_eval, .dtor = inter_free};
  sp_func_t* intersection = sp_func_new(inter_str, inter, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (has_grad)
  {
    char inter_grad_str[1024];
    sprintf(inter_grad_str, "Intersection gradient"); // FIXME: Yadda
    sp_vtable vtable_g = {.eval = inter_eval_gradient}; // Notice no dtor.
    sp_func_t* inter_grad = sp_func_new(inter_grad_str, inter, vtable_g, SP_INHOMOGENEOUS, 3);
    sp_func_register_deriv(intersection, 1, inter_grad);
  }

  return intersection;
}

sp_func_t* intersection_new2(sp_func_t* surface1, sp_func_t* surface2)
{
  sp_func_t* surfs[2];
  surfs[0] = surface1;
  surfs[1] = surface2;
  return intersection_new(surfs, 2);
}

