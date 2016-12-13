// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sp_func.h"

struct sp_func_t 
{
  char* name;
  void* context;
  sp_func_vtable vtable;
  int num_comp;
  bool homogeneous;
};

static void sp_func_free(void* ctx)
{
  sp_func_t* func = ctx;
  if ((func->vtable.dtor != NULL) && (func->context != NULL))
    func->vtable.dtor(func->context);
  func->context = NULL;
  polymec_free(func->name);
}

sp_func_t* sp_func_new(const char* name, void* context, sp_func_vtable vtable,
                       sp_func_homogeneity_t homogeneity,
                       int num_comp)
{
  ASSERT(context != NULL);
  ASSERT(vtable.eval != NULL);
  ASSERT(num_comp > 0);
  sp_func_t* f = polymec_gc_malloc(sizeof(sp_func_t), sp_func_free);
  f->name = string_dup(name);
  f->context = context;
  f->vtable = vtable;
  f->homogeneous = (homogeneity == SP_FUNC_HOMOGENEOUS);
  f->num_comp = num_comp;
  return f;
}

sp_func_t* sp_func_from_func(const char* name, sp_eval_func func, 
                             sp_func_homogeneity_t homogeneity,
                             int num_comp)
{
  ASSERT(func != NULL);
  ASSERT(num_comp > 0);
  sp_func_t* f = polymec_gc_malloc(sizeof(sp_func_t), sp_func_free);
  f->name = string_dup(name);
  f->context = NULL;
  f->vtable.eval = func;
  f->homogeneous = (homogeneity == SP_FUNC_HOMOGENEOUS);
  f->num_comp = num_comp;
  return f;
}

const char* sp_func_name(sp_func_t* func)
{
  return (const char*)func->name;
}

void sp_func_rename(sp_func_t* func, const char* new_name)
{
  polymec_free(func->name);
  func->name = string_dup(new_name);
}

bool sp_func_is_homogeneous(sp_func_t* func)
{
  return func->homogeneous;
}

int sp_func_num_comp(sp_func_t* func)
{
  return func->num_comp;
}

void* sp_func_context(sp_func_t* func)
{
  return func->context;
}

void sp_func_eval(sp_func_t* func, point_t* x, real_t* result)
{
  func->vtable.eval(func->context, x, result);
}

// Constant spatial function implementation.

typedef struct
{
  int num_comp;
  real_t *comp;
} const_sp_func_t;

static void constant_eval(void* ctx, point_t* x, real_t* res)
{
  const_sp_func_t* f = ctx;
  for (int i = 0; i < f->num_comp; ++i)
    res[i] = f->comp[i];
}

static void constant_dtor(void* ctx)
{
  const_sp_func_t* f = ctx;
  polymec_free(f->comp);
  polymec_free(f);
}

sp_func_t* constant_sp_func_new(real_t components[], int num_components)
{
  sp_func_vtable vtable = {.eval = constant_eval, .dtor = constant_dtor};
  char name[1025];
  snprintf(name, 1024, "constant spatial function (");
  for (int i = 0; i < num_components; ++i)
  {
    char comp_str[22];
    if (i == (num_components-1))
      snprintf(comp_str, 20, "%g)", components[i]);
    else
      snprintf(comp_str, 20, "%g, ", components[i]);
    strncat(name, comp_str, 1024);
  }
  const_sp_func_t* f = polymec_malloc(sizeof(const_sp_func_t));
  f->num_comp = num_components;
  f->comp = polymec_malloc(num_components*sizeof(real_t));
  for (int i = 0; i < num_components; ++i)
    f->comp[i] = components[i];
  return sp_func_new(name, (void*)f, vtable, SP_FUNC_HOMOGENEOUS, num_components);
}

