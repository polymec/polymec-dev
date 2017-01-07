// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/sd_func.h"

struct sd_func_t 
{
  char* name;
  void* context;
  sd_func_vtable vtable;
};

static void sd_func_free(void* ctx)
{
  sd_func_t* func = ctx;
  if ((func->vtable.dtor != NULL) && (func->context != NULL))
    func->vtable.dtor(func->context);
  func->context = NULL;
  polymec_free(func->name);
}

sd_func_t* sd_func_new(const char* name, void* context, sd_func_vtable vtable)
{
  ASSERT(vtable.value != NULL);
  ASSERT(vtable.eval_grad != NULL);
  sd_func_t* f = polymec_gc_malloc(sizeof(sd_func_t), sd_func_free);
  f->name = string_dup(name);
  f->context = context;
  f->vtable = vtable;
  return f;
}

typedef struct 
{
  sp_func_t* f;
  sp_func_t* df;
} spsd_t;

static real_t spsd_value(void* context, point_t* x)
{
  spsd_t* spsd = context;
  real_t val;
  sp_func_eval(spsd->f, x, &val);
  return val;
}

static void spsd_eval_grad(void* context, point_t* x, vector_t* grad)
{
  spsd_t* spsd = context;
  sp_func_eval(spsd->df, x, (real_t*)grad);
}

static void spsd_dtor(void* context)
{
  spsd_t* spsd = context;
  spsd->f = NULL;
  spsd->df = NULL;
  polymec_free(spsd);
}

sd_func_t* sd_func_from_sp_funcs(const char* name, 
                                 sp_func_t* distance, 
                                 sp_func_t* gradient)
{
  ASSERT(distance != NULL);
  ASSERT(gradient != NULL);
  ASSERT(sp_func_num_comp(distance) == 1);
  ASSERT(sp_func_num_comp(gradient) == 3);
  spsd_t* spsd = polymec_malloc(sizeof(spsd_t));
  spsd->f = distance;
  spsd->df = gradient;
  sd_func_vtable vtable = {.value = spsd_value,
                           .eval_grad = spsd_eval_grad,
                           .dtor = spsd_dtor};
  return sd_func_new(name, spsd, vtable);
}

const char* sd_func_name(sd_func_t* func)
{
  return (const char*)func->name;
}

void sd_func_rename(sd_func_t* func, const char* new_name)
{
  polymec_free(func->name);
  func->name = string_dup(new_name);
}

void* sd_func_context(sd_func_t* func)
{
  return func->context;
}

real_t sd_func_value(sd_func_t* func, point_t* x)
{
  return func->vtable.value(func->context, x);
}

void sd_func_eval_grad(sd_func_t* func, point_t* x, vector_t* grad)
{
  func->vtable.eval_grad(func->context, x, grad);
}

void sd_func_project(sd_func_t* func, point_t* x, point_t* proj_x)
{
  real_t D = func->vtable.value(func->context, x);
  vector_t grad;
  func->vtable.eval_grad(func->context, x, &grad);
  real_t G = vector_mag(&grad);
  if (G > 0.0)
  {
    proj_x->x = x->x - D * grad.x / G;
    proj_x->y = x->y - D * grad.y / G;
    proj_x->z = x->z - D * grad.z / G;
  }
  else
    *proj_x = *x;
}
