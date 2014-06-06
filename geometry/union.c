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

#include "geometry/union.h"

typedef struct
{
  sp_func_t** funcs;
  int num_funcs;
} un_t;

// Destructor.
static void un_free(void* ctx)
{
  un_t* un = ctx;
  polymec_free(un->funcs);
  polymec_free(un);
}

static void un_eval(void* ctx, point_t* x, real_t* result)
{
  un_t* un = ctx;
  real_t minval = FLT_MAX;
  for (int i = 0; i < un->num_funcs; ++i)
  {
    real_t ival;
    sp_func_eval(un->funcs[i], x, &ival);
    minval = MIN(minval, ival);
  }
  result[0] = minval;
}

static void un_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  un_t* un = ctx;
  real_t minval = FLT_MAX;
  int index = -1;
  for (int i = 0; i < un->num_funcs; ++i)
  {
    real_t ival;
    sp_func_eval(un->funcs[i], x, &ival);
    if (ival < minval)
    {
      minval = ival;
      index = i;
    }
  }
  sp_func_eval(un->funcs[index], x, result);
}

sp_func_t* union_new(sp_func_t** surfaces, int num_surfaces)
{
  ASSERT(surfaces != NULL);
  ASSERT(num_surfaces > 1);

  un_t* un = polymec_malloc(sizeof(un_t));
  un->funcs = polymec_malloc(sizeof(sp_func_t*) * num_surfaces);
  un->num_funcs = num_surfaces;
  bool has_grad = true;
  for (int i = 0; i < num_surfaces; ++i)
  {
    ASSERT(sp_func_num_comp(surfaces[i]) == 1);
    un->funcs[i] = surfaces[i];
    if (!sp_func_has_deriv(surfaces[i], 1))
      has_grad = false;
  }

  char un_str[num_surfaces*1024+1];
  sprintf(un_str, "union (");
  for (int i = 0; i < num_surfaces; ++i)
  {
    char surf_str[1024];
    if (i == (num_surfaces - 1))
      snprintf(surf_str, 1024, "%s)", sp_func_name(surfaces[i]));
    else
      snprintf(surf_str, 1024, "%s, ", sp_func_name(surfaces[i]));
    strncat(un_str, surf_str, num_surfaces*1024);
  }
  sp_vtable vtable = {.eval = un_eval, .dtor = un_free};
  sp_func_t* union_func = sp_func_new(un_str, un, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (has_grad)
  {
    char un_grad_str[num_surfaces*1024+1];
    snprintf(un_grad_str, num_surfaces*1024, "grad %s", un_str);
    sp_vtable vtable_g = {.eval = un_eval_gradient}; // Notice no dtor.
    sp_func_t* un_grad = sp_func_new(un_grad_str, un, vtable_g, SP_INHOMOGENEOUS, 3);
    sp_func_register_deriv(union_func, 1, un_grad);
  }

  return union_func;
}

sp_func_t* union_new2(sp_func_t* surface1, sp_func_t* surface2)
{
  sp_func_t* surfs[2];
  surfs[0] = surface1;
  surfs[1] = surface2;
  return union_new(surfs, 2);
}

