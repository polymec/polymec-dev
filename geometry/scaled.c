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

#include "geometry/scaled.h"

typedef struct
{
  sp_func_t* func;
  double scale_factor;
} sc_t;

static void sc_eval(void* ctx, point_t* x, double* result)
{
  sc_t* sc = ctx;
  // Scale x.
  double A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sp_func_eval(sc->func, &y, result);
}

static void sc_eval_gradient(void* ctx, point_t* x, double* result)
{
  sc_t* sc = ctx;
  ASSERT(sp_func_has_deriv(sc->func, 1));

  // Scale x.
  double A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sp_func_eval_deriv(sc->func, 1, &y, result);
}

sp_func_t* scaled_new(sp_func_t* func, double scale_factor)
{
  ASSERT(func != NULL);
  ASSERT(sp_func_num_comp(func) == 1);
  ASSERT(scale_factor > 0.0);

  sc_t* sc = malloc(sizeof(sc_t));
  sc->func = func;
  sc->scale_factor = scale_factor;
  char sc_str[1024];
  sprintf(sc_str, "scaled"); // FIXME: Not very helpful.
  sp_vtable vtable = {.eval = sc_eval};
  sp_func_t* sc_func = sp_func_new(sc_str, sc, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (sp_func_has_deriv(func, 1))
  {
    char sc_grad_str[1024];
    sprintf(sc_grad_str, "scaled gradient"); // FIXME: Yadda
    sp_vtable vtable_g = {.eval = sc_eval_gradient}; // Notice no dtor.
    sp_func_t* sc_grad = sp_func_new(sc_grad_str, sc, vtable_g, SP_INHOMOGENEOUS, 3);
    sp_func_register_deriv(sc_func, 1, sc_grad);
  }

  return sc_func;
}

