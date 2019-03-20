// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
  ASSERT(vtable.eval_grad != NULL);
  sdt_func_t* f = polymec_refcounted_malloc(sizeof(sdt_func_t), sdt_func_free);
  f->name = string_dup(name);
  f->context = context;
  f->vtable = vtable;
  return f;
}

typedef struct
{
  st_func_t* f;
  st_func_t* df;
} stsd_t;

static real_t stsd_value(void* context, point_t* x, real_t t)
{
  stsd_t* stsd = context;
  real_t val;
  st_func_eval(stsd->f, x, t, &val);
  return val;
}

static void stsd_eval_n(void* context, point_t* xs, size_t n, real_t t, real_t* vals)
{
  stsd_t* stsd = context;
  st_func_eval_n(stsd->f, xs, n, t, vals);
}


static void stsd_eval_grad(void* context, point_t* x, real_t t, vector_t* grad)
{
  stsd_t* stsd = context;
  st_func_eval(stsd->df, x, t, (real_t*)grad);
}

static void stsd_eval_n_grad(void* context, point_t* xs, size_t n, real_t t, vector_t* grads)
{
  stsd_t* stsd = context;
  st_func_eval_n(stsd->df, xs, n, t, (real_t*)grads);
}

static void stsd_dtor(void* context)
{
  stsd_t* stsd = context;
  stsd->f = NULL;
  stsd->df = NULL;
  polymec_free(stsd);
}

sdt_func_t* sdt_func_from_st_funcs(const char* name,
                                   st_func_t* distance,
                                   st_func_t* gradient)
{
  ASSERT(distance != NULL);
  ASSERT(gradient != NULL);
  ASSERT(st_func_num_comp(distance) == 1);
  ASSERT(st_func_num_comp(gradient) == 3);
  stsd_t* stsd = polymec_malloc(sizeof(stsd_t));
  stsd->f = distance;
  stsd->df = gradient;
  sdt_func_vtable vtable = {.value = stsd_value,
                            .eval_n = stsd_eval_n,
                            .eval_grad = stsd_eval_grad,
                            .eval_n_grad = stsd_eval_n_grad,
                            .dtor = stsd_dtor};
  return sdt_func_new(name, stsd, vtable);
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

void sdt_func_eval_n(sdt_func_t* func, point_t* xs, size_t n, real_t t, real_t* vals)
{
  if (func->vtable.eval_n != NULL)
    func->vtable.eval_n(func->context, xs, n, t, vals);
  else
  {
    for (size_t i = 0; i < n; ++i)
      vals[i] = func->vtable.value(func->context, &xs[i], t);
  }
}

void sdt_func_eval_grad(sdt_func_t* func, point_t* x, real_t t, vector_t* grad)
{
  func->vtable.eval_grad(func->context, x, t, grad);
}

void sdt_func_eval_n_grad(sdt_func_t* func, point_t* xs, size_t n, real_t t, vector_t* grads)
{
  if (func->vtable.eval_n_grad != NULL)
    func->vtable.eval_n_grad(func->context, xs, n, t, grads);
  else
  {
    for (size_t i = 0; i < n; ++i)
      func->vtable.eval_grad(func->context, &xs[i], t, &grads[i]);
  }
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

void sdt_func_project_n(sdt_func_t* func, point_t* xs, size_t n, real_t t, point_t* proj_xs)
{
  real_t* Ds = polymec_malloc(sizeof(real_t) * n);
  sdt_func_eval_n(func, xs, n, t, Ds);
  vector_t* grads = polymec_malloc(sizeof(vector_t) * n);
  sdt_func_eval_n_grad(func->context, xs, n, t, grads);
  for (size_t i = 0; i < n; ++i)
  {
    real_t D = Ds[i];
    real_t G = vector_mag(&grads[i]);
    if (G > 0.0)
    {
      proj_xs[i].x = xs[i].x - D * grads[i].x / G;
      proj_xs[i].y = xs[i].y - D * grads[i].y / G;
      proj_xs[i].z = xs[i].z - D * grads[i].z / G;
    }
    else
      proj_xs[i] = xs[i];
  }
  polymec_free(grads);
  polymec_free(Ds);
}
