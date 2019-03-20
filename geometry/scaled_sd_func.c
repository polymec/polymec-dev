// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/scaled_sd_func.h"

typedef struct
{
  sd_func_t* func;
  real_t scale_factor;
} sc_t;

static real_t sc_value(void* ctx, point_t* x)
{
  sc_t* sc = ctx;
  real_t A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  return sd_func_value(sc->func, &y);
}

static void sc_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  sc_t* sc = ctx;
  real_t A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sd_func_eval_grad(sc->func, &y, grad);
}

sd_func_t* scaled_sd_func_new(sd_func_t* func, real_t scale_factor)
{
  ASSERT(func != NULL);
  ASSERT(scale_factor > 0.0);

  sc_t* sc = polymec_malloc(sizeof(sc_t));
  sc->func = func;
  sc->scale_factor = scale_factor;
  char sc_str[1024];
  snprintf(sc_str, 1024, "scaled(%s)", sd_func_name(func));
  sd_func_vtable vtable = {.value = sc_value,
                           .eval_grad = sc_eval_grad,
                           .dtor = polymec_free};
  return sd_func_new(sc_str, sc, vtable);
}

