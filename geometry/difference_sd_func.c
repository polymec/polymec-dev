// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/difference_sd_func.h"

typedef struct
{
  sd_func_t *s1, *s2;
} diff_t;

static real_t diff_value(void* ctx, point_t* x)
{
  diff_t* diff = ctx;
  real_t v1 = sd_func_value(diff->s1, x);
  real_t v2 = sd_func_value(diff->s2, x);
  return MAX(v1, -v2);
}

static void diff_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  diff_t* diff = ctx;
  real_t v1 = sd_func_value(diff->s1, x);
  real_t v2 = sd_func_value(diff->s2, x);
  if (v1 > v2)
    sd_func_eval_grad(diff->s1, x, grad);
  else
    sd_func_eval_grad(diff->s2, x, grad);
}

sd_func_t* difference_sd_func_new(sd_func_t* surface1, sd_func_t* surface2)
{
  ASSERT(surface1 != NULL);
  ASSERT(surface2 != NULL);

  diff_t* diff = polymec_malloc(sizeof(diff_t));
  diff->s1 = surface1;
  diff->s2 = surface2;

  char diff_str[4096];
  sprintf(diff_str, "Difference (%s, %s)", sd_func_name(surface1), sd_func_name(surface2));
  sd_func_vtable vtable = {.value = diff_value,
                           .eval_grad = diff_eval_grad,
                           .dtor = polymec_free};
  return sd_func_new(diff_str, diff, vtable);
}

