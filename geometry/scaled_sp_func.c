// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/scaled_sp_func.h"

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

sp_func_t* scaled_sp_func_new(sp_func_t* func, real_t scale_factor)
{
  ASSERT(func != NULL);
  ASSERT(sp_func_num_comp(func) == 1);
  ASSERT(scale_factor > 0.0);

  sc_t* sc = polymec_malloc(sizeof(sc_t));
  sc->func = func;
  sc->scale_factor = scale_factor;
  char sc_str[1024];
  snprintf(sc_str, 1024, "scaled(%s)", sp_func_name(func)); 
  sp_func_vtable vtable = {.eval = sc_eval, .dtor = polymec_free};
  sp_func_t* sc_func = sp_func_new(sc_str, sc, vtable, SP_FUNC_HETEROGENEOUS, 1);

  // Register the gradient function if we have it.
  if (sp_func_has_deriv(func, 1))
  {
    char sc_grad_str[1024];
    snprintf(sc_grad_str, 1024, "grad %s", sc_str);
    sp_func_vtable vtable_g = {.eval = sc_eval_gradient}; // Notice no dtor.
    sp_func_t* sc_grad = sp_func_new(sc_grad_str, sc, vtable_g, SP_FUNC_HETEROGENEOUS, 3);
    sp_func_register_deriv(sc_func, 1, sc_grad);
  }

  return sc_func;
}

