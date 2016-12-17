// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/sdt_func.h"

struct sdt_func_t 
{
  char* name;
  void* context;
  sdt_func_vtable vtable;
};

static void sdt_func_free(void* ctx)
{
  sdt_func_t* func = ctx;
  if ((func->vtable.dtor != NULL) && (func->context != NULL))
    func->vtable.dtor(func->context);
  func->context = NULL;
  polymec_free(func->name);
}

sdt_func_t* sdt_func_new(const char* name, void* context, sdt_func_vtable vtable)
{
  ASSERT(vtable.value != NULL);
  ASSERT(vtable.eval_grad == NULL);
  sdt_func_t* f = polymec_gc_malloc(sizeof(sdt_func_t), sdt_func_free);
  f->name = string_dup(name);
  f->context = context;
  f->vtable = vtable;
  return f;
}

typedef struct 
{
  st_func_t* f;
  st_func_t* df;
} spsd_t;

static real_t spsd_value(void* context, point_t* x, real_t t)
{
  spsd_t* spsd = context;
  real_t val;
  st_func_eval(spsd->f, x, t, &val);
  return val;
}

static void spsd_eval_grad(void* context, point_t* x, real_t t, vector_t* grad)
{
  spsd_t* spsd = context;
  st_func_eval(spsd->df, x, t, (real_t*)grad);
}

static void spsd_dtor(void* context)
{
  spsd_t* spsd = context;
  spsd->f = NULL;
  spsd->df = NULL;
  polymec_free(spsd);
}

sdt_func_t* sdt_func_from_st_funcs(const char* name, 
                                   st_func_t* distance, 
                                   st_func_t* gradient)
{
  ASSERT(distance != NULL);
  ASSERT(gradient != NULL);
  ASSERT(st_func_num_comp(distance) == 1);
  ASSERT(st_func_num_comp(gradient) == 3);
  spsd_t* spsd = polymec_malloc(sizeof(spsd_t));
  spsd->f = distance;
  spsd->df = gradient;
  sdt_func_vtable vtable = {.value = spsd_value,
                            .eval_grad = spsd_eval_grad,
                            .dtor = spsd_dtor};
  return sdt_func_new(name, spsd, vtable);
}

const char* sdt_func_name(sdt_func_t* func)
{
  return (const char*)func->name;
}

void sdt_func_rename(sdt_func_t* func, const char* new_name)
{
  polymec_free(func->name);
  func->name = string_dup(new_name);
}

void* sdt_func_context(sdt_func_t* func)
{
  return func->context;
}

real_t sdt_func_value(sdt_func_t* func, point_t* x, real_t t)
{
  return func->vtable.value(func->context, x, t);
}

void sdt_func_eval_grad(sdt_func_t* func, point_t* x, real_t t, vector_t* grad)
{
  func->vtable.eval_grad(func->context, x, t, grad);
}

void sdt_func_project(sdt_func_t* func, point_t* x, real_t t, point_t* proj_x)
{
  real_t D = func->vtable.value(func->context, x, t);
  vector_t grad;
  func->vtable.eval_grad(func->context, x, t, &grad);
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
