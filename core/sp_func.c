// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "core/sp_func.h"

struct sp_func_t 
{
  char* name;
  void* context;
  sp_vtable vtable;
  int num_comp;
  bool homogeneous;

  // Functions for computing derivatives (up to 4).
  sp_func_t* derivs[4];
};

static void sp_func_free(void* ctx, void* dummy)
{
  sp_func_t* func = ctx;
  if ((func->vtable.dtor != NULL) && (func->context != NULL))
    func->vtable.dtor(func->context);
  func->context = NULL;
  polymec_free(func->name);
}

sp_func_t* sp_func_new(const char* name, void* context, sp_vtable vtable,
                       sp_func_homogeneity_t homogeneity,
                       int num_comp)
{
  ASSERT(context != NULL);
  ASSERT(vtable.eval != NULL);
  ASSERT((vtable.eval_deriv == NULL) || (vtable.has_deriv != NULL));
  ASSERT(num_comp > 0);
  sp_func_t* f = GC_MALLOC(sizeof(sp_func_t));
  f->name = string_dup(name);
  f->context = context;
  f->vtable = vtable;
  f->homogeneous = (homogeneity == SP_HOMOGENEOUS);
  f->num_comp = num_comp;
  memset(f->derivs, 0, sizeof(sp_func_t*)*4);
  GC_register_finalizer(f, sp_func_free, f, NULL, NULL);
  return f;
}

sp_func_t* sp_func_from_func(const char* name, sp_eval_func func, 
                             sp_func_homogeneity_t homogeneity,
                             int num_comp)
{
  ASSERT(func != NULL);
  ASSERT(num_comp > 0);
  sp_func_t* f = GC_MALLOC(sizeof(sp_func_t));
  f->name = string_dup(name);
  f->context = NULL;
  f->vtable.eval = func;
  f->homogeneous = (homogeneity == SP_HOMOGENEOUS);
  f->num_comp = num_comp;
  memset(f->derivs, 0, sizeof(sp_func_t*)*4);
  GC_register_finalizer(f, &sp_func_free, f, NULL, NULL);
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

void sp_func_register_deriv(sp_func_t* func, int n, sp_func_t* nth_deriv)
{
  ASSERT(n > 0);
  ASSERT(n <= 4);
  ASSERT(nth_deriv != NULL);
  ASSERT(sp_func_num_comp(nth_deriv) == (func->num_comp * (int)pow(3, n))); 
  func->derivs[n-1] = nth_deriv;
}

bool sp_func_has_deriv(sp_func_t* func, int n)
{
  ASSERT(n > 0);
  if (n > 4) return false; // FIXME
  if (func->vtable.eval_deriv != NULL)
    return func->vtable.has_deriv(func->context, n);
  else
    return (func->derivs[n-1] != NULL);
}

// Evaluates the nth derivative of this function, placing the result in result.
void sp_func_eval_deriv(sp_func_t* func, int n, point_t* x, real_t* result)
{
  if (func->vtable.eval_deriv != NULL)
    func->vtable.eval_deriv(func->context, n, x, result);
  else
    sp_func_eval(func->derivs[n-1], x, result);
}

#if 0
void sp_func_grad_centered_diff(sp_func_t* func, point_t* x0, vector_t* dx, vector_t* gradient)
{
  point_t x1 = {x0->x-dx->x, x0->y, x0->z}, 
          x2 = {x0->x+dx->x, x0->y, x0->z};
  real_t f1, f2;
  sp_func_eval(func, &x1, &f1);
  sp_func_eval(func, &x2, &f2);
  gradient->x = (f2-f1) / dx->x;
  x1.x = x0->x, x2.x = x0->x;
  x1.y = x0->y-dx->y, x2.y = x0->y+dx->y;
  sp_func_eval(func, &x1, &f1);
  sp_func_eval(func, &x2, &f2);
  gradient->y = (f2-f1) / dx->y;
  x1.y = x0->y, x2.y = x0->y;
  x1.z = x0->z-dx->z, x2.z = x0->z+dx->z;
  sp_func_eval(func, &x1, &f1);
  sp_func_eval(func, &x2, &f2);
  gradient->z = (f2-f1) / dx->z;
}

void sp_func_grad_richardson(sp_func_t* func, point_t* x0, vector_t* dx, vector_t* gradient)
{
  // Form the low-res finite difference approximation of the gradient.
  vector_t Gl;
  sp_func_grad_centered_diff(func, x0, dx, &Gl);

  // Form the hi-res finite difference approximation of the gradient.
  vector_t Gh;
  vector_t dx1 = {0.5*dx->x, 0.5*dx->y, 0.5*dx->z};
  sp_func_grad_centered_diff(func, x0, &dx1, &Gh);

  // Form the Richardson extrapolation.
  static int N = 2;
  real_t twoN = pow(2.0, N);
  gradient->x = (twoN * Gh.x - Gl.x) / (twoN - 1.0);
  gradient->y = (twoN * Gh.y - Gl.y) / (twoN - 1.0);
  gradient->z = (twoN * Gh.z - Gl.z) / (twoN - 1.0);
}
#endif

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

static sp_func_t* create_constant_sp_func(real_t components[], int num_components)
{
  sp_vtable vtable = {.eval = constant_eval, .dtor = constant_dtor};
  char name[1024];
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
  return sp_func_new(name, (void*)f, vtable, SP_HOMOGENEOUS, num_components);
}

sp_func_t* constant_sp_func_new(real_t components[], int num_components)
{
  sp_func_t* func = create_constant_sp_func(components, num_components);

  // Just to be complete, we register the 3-component zero function as this 
  // function's gradient.
  real_t zeros[3*num_components];
  memset(zeros, 0, 3*num_components*sizeof(real_t));
  sp_func_t* zero = create_constant_sp_func(zeros, 3*num_components);
  sp_func_register_deriv(func, 1, zero);

  return func;
}


