// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/intersection_sp_func.h"

typedef struct
{
  sp_func_t** funcs;
  int num_funcs;
} inter_t;

// Destructor.
static void inter_free(void* ctx)
{
  inter_t* inter = ctx;
  polymec_free(inter->funcs);
  polymec_free(inter);
}

static void inter_eval(void* ctx, point_t* x, real_t* result)
{
  inter_t* inter = ctx;
  real_t maxval = -FLT_MAX;
  for (int i = 0; i < inter->num_funcs; ++i)
  {
    real_t ival;
    sp_func_eval(inter->funcs[i], x, &ival);
    maxval = MAX(maxval, ival);
  }
  result[0] = maxval;
}

static void inter_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  inter_t* inter = ctx;
  real_t maxval = -FLT_MAX;
  int index = -1;
  for (int i = 0; i < inter->num_funcs; ++i)
  {
    real_t ival;
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

sp_func_t* intersection_sp_func_new(sp_func_t** surfaces, int num_surfaces)
{
  ASSERT(surfaces != NULL);
  ASSERT(num_surfaces > 1);

  inter_t* inter = polymec_malloc(sizeof(inter_t));
  inter->funcs = polymec_malloc(sizeof(sp_func_t*) * num_surfaces);
  inter->num_funcs = num_surfaces;
  bool has_grad = true;
  for (int i = 0; i < num_surfaces; ++i)
  {
    ASSERT(sp_func_num_comp(surfaces[i]) == 1);
    inter->funcs[i] = surfaces[i];
    if (!sp_func_has_deriv(surfaces[i], 1))
      has_grad = false;
  }

  char inter_str[num_surfaces*1024+1];
  sprintf(inter_str, "union (");
  for (int i = 0; i < num_surfaces; ++i)
  {
    char surf_str[1024];
    if (i == (num_surfaces - 1))
      snprintf(surf_str, 1024, "%s)", sp_func_name(surfaces[i]));
    else
      snprintf(surf_str, 1024, "%s, ", sp_func_name(surfaces[i]));
    strncat(inter_str, surf_str, num_surfaces*1024);
  }
  sp_func_vtable vtable = {.eval = inter_eval, .dtor = inter_free};
  sp_func_t* intersection = sp_func_new(inter_str, inter, vtable, SP_FUNC_HETEROGENEOUS, 1);

  // Register the gradient function if we have it.
  if (has_grad)
  {
    char inter_grad_str[num_surfaces*1024];
    snprintf(inter_grad_str, num_surfaces*1024, "grad %s", inter_str);
    sp_func_vtable vtable_g = {.eval = inter_eval_gradient}; // Notice no dtor.
    sp_func_t* inter_grad = sp_func_new(inter_grad_str, inter, vtable_g, SP_FUNC_HETEROGENEOUS, 3);
    sp_func_register_deriv(intersection, 1, inter_grad);
  }

  return intersection;
}

sp_func_t* intersection_sp_func_new2(sp_func_t* surface1, sp_func_t* surface2)
{
  sp_func_t* surfs[2];
  surfs[0] = surface1;
  surfs[1] = surface2;
  return intersection_sp_func_new(surfs, 2);
}

