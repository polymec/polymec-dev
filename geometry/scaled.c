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

#include "geometry/scaled.h"

typedef struct
{
  sp_func_t* func;
  real_t scale_factor;
} sc_t;

static void sc_eval(void* ctx, point_t* x, real_t* result)
{
  sc_t* sc = ctx;
  // Scale x.
  real_t A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sp_func_eval(sc->func, &y, result);
}

static void sc_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  sc_t* sc = ctx;
  ASSERT(sp_func_has_deriv(sc->func, 1));

  // Scale x.
  real_t A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sp_func_eval_deriv(sc->func, 1, &y, result);
}

sp_func_t* scaled_new(sp_func_t* func, real_t scale_factor)
{
  ASSERT(func != NULL);
  ASSERT(sp_func_num_comp(func) == 1);
  ASSERT(scale_factor > 0.0);

  sc_t* sc = polymec_malloc(sizeof(sc_t));
  sc->func = func;
  sc->scale_factor = scale_factor;
  char sc_str[1024];
  sprintf(sc_str, "scaled"); // FIXME: Not very helpful.
  sp_vtable vtable = {.eval = sc_eval, .dtor = polymec_free};
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

