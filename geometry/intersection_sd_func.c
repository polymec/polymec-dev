// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/intersection_sd_func.h"

typedef struct
{
  sd_func_t** funcs;
  int num_funcs;
} inter_t;

// Destructor.
static void inter_free(void* ctx)
{
  inter_t* inter = ctx;
  polymec_free(inter->funcs);
  polymec_free(inter);
}

static real_t inter_value(void* ctx, point_t* x)
{
  inter_t* inter = ctx;
  real_t maxval = -REAL_MAX;
  for (int i = 0; i < inter->num_funcs; ++i)
    maxval = MAX(maxval, sd_func_value(inter->funcs[i], x));
  return maxval;
}

static void inter_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  inter_t* inter = ctx;
  real_t maxval = -REAL_MAX;
  int index = -1;
  for (int i = 0; i < inter->num_funcs; ++i)
  {
    real_t ival = sd_func_value(inter->funcs[i], x);
    if (ival > maxval)
    {
      maxval = ival;
      index = i;
    }
  }
  sd_func_eval_grad(inter->funcs[index], x, grad);
}

sd_func_t* intersection_sd_func_new(sd_func_t** surfaces, int num_surfaces)
{
  ASSERT(surfaces != NULL);
  ASSERT(num_surfaces > 1);

  inter_t* inter = polymec_malloc(sizeof(inter_t));
  inter->funcs = polymec_malloc(sizeof(sd_func_t*) * num_surfaces);
  inter->num_funcs = num_surfaces;
  for (int i = 0; i < num_surfaces; ++i)
    inter->funcs[i] = surfaces[i];

  char inter_str[num_surfaces*1024+1];
  sprintf(inter_str, "union (");
  for (int i = 0; i < num_surfaces; ++i)
  {
    char surf_str[1024];
    if (i == (num_surfaces - 1))
      snprintf(surf_str, 1024, "%s)", sd_func_name(surfaces[i]));
    else
      snprintf(surf_str, 1024, "%s, ", sd_func_name(surfaces[i]));
    strncat(inter_str, surf_str, num_surfaces*1024);
  }
  sd_func_vtable vtable = {.value = inter_value,
                           .eval_grad = inter_eval_grad,
                           .dtor = inter_free};
  return sd_func_new(inter_str, inter, vtable);
}

sd_func_t* intersection_sd_func_new2(sd_func_t* surface1, sd_func_t* surface2)
{
  sd_func_t* surfs[2];
  surfs[0] = surface1;
  surfs[1] = surface2;
  return intersection_sd_func_new(surfs, 2);
}

