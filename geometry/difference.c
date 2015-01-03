// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/difference.h"

typedef struct
{
  sp_func_t *s1, *s2;
} diff_t;

static void diff_eval(void* ctx, point_t* x, real_t* result)
{
  diff_t* diff = ctx;
  real_t v1, v2;
  sp_func_eval(diff->s1, x, &v1);
  sp_func_eval(diff->s2, x, &v2);
  result[0] = MAX(v1, -v2);
}

static void diff_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  diff_t* diff = ctx;
  real_t v1[3], v2[3];
  sp_func_eval(diff->s1, x, v1);
  real_t v1mag2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
  sp_func_eval(diff->s2, x, v2);
  real_t v2mag2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];

  if (v1mag2 > v2mag2)
  {
    result[0] = v1[0];
    result[1] = v1[1];
    result[2] = v1[2];
  }
  else
  {
    result[0] = v2[0];
    result[1] = v2[1];
    result[2] = v2[2];
  }
}

sp_func_t* difference_new(sp_func_t* surface1, sp_func_t* surface2)
{
  ASSERT(surface1 != NULL);
  ASSERT(sp_func_num_comp(surface1) == 1);
  ASSERT(surface2 != NULL);
  ASSERT(sp_func_num_comp(surface2) == 1);

  diff_t* diff = polymec_malloc(sizeof(diff_t));
  diff->s1 = surface1;
  diff->s2 = surface2;

  char diff_str[4096];
  sprintf(diff_str, "Difference (%s, %s)", sp_func_name(surface1), sp_func_name(surface2));
  sp_vtable vtable = {.eval = diff_eval, .dtor = polymec_free};
  sp_func_t* difference = sp_func_new(diff_str, diff, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (sp_func_has_deriv(surface1, 1) && sp_func_has_deriv(surface2, 1))
  {
    char diff_grad_str[4096];
    sprintf(diff_grad_str, "grad %s", diff_str);
    sp_vtable vtable_g = {.eval = diff_eval_gradient}; // Notice no dtor.
    sp_func_t* diff_grad = sp_func_new(diff_grad_str, diff, vtable_g, SP_INHOMOGENEOUS, 3);
    sp_func_register_deriv(difference, 1, diff_grad);
  }

  return difference;
}

