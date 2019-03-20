// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/union_sd_func.h"

typedef struct
{
  sd_func_t** funcs;
  int num_funcs;
} un_t;

// Destructor.
static void un_free(void* ctx)
{
  un_t* un = ctx;
  polymec_free(un->funcs);
  polymec_free(un);
}

static real_t un_value(void* ctx, point_t* x)
{
  un_t* un = ctx;
  real_t minval = REAL_MAX;
  for (int i = 0; i < un->num_funcs; ++i)
    minval = MIN(minval, sd_func_value(un->funcs[i], x));
  return minval;
}

static void un_eval_grad(void* ctx, point_t* x, vector_t* grad)
{
  un_t* un = ctx;
  real_t minval = REAL_MAX;
  int index = -1;
  for (int i = 0; i < un->num_funcs; ++i)
  {
    real_t ival = sd_func_value(un->funcs[i], x);
    if (ival < minval)
    {
      minval = ival;
      index = i;
    }
  }
  sd_func_eval_grad(un->funcs[index], x, grad);
}

sd_func_t* union_sd_func_new(sd_func_t** surfaces, int num_surfaces)
{
  ASSERT(surfaces != NULL);
  ASSERT(num_surfaces > 1);

  un_t* un = polymec_malloc(sizeof(un_t));
  un->funcs = polymec_malloc(sizeof(sd_func_t*) * num_surfaces);
  un->num_funcs = num_surfaces;
  for (int i = 0; i < num_surfaces; ++i)
    un->funcs[i] = surfaces[i];

  char un_str[num_surfaces*1024+1];
  sprintf(un_str, "union (");
  for (int i = 0; i < num_surfaces; ++i)
  {
    char surf_str[1024];
    if (i == (num_surfaces - 1))
      snprintf(surf_str, 1024, "%s)", sd_func_name(surfaces[i]));
    else
      snprintf(surf_str, 1024, "%s, ", sd_func_name(surfaces[i]));
    strncat(un_str, surf_str, num_surfaces*1024);
  }
  sd_func_vtable vtable = {.value = un_value,
                           .eval_grad = un_eval_grad,
                           .dtor = un_free};
  return sd_func_new(un_str, un, vtable);
}

sd_func_t* union_sd_func_new2(sd_func_t* surface1, sd_func_t* surface2)
{
  sd_func_t* surfs[2];
  surfs[0] = surface1;
  surfs[1] = surface2;
  return union_sd_func_new(surfs, 2);
}

