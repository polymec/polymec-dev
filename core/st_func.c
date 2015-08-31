// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "core/st_func.h"

struct st_func_t 
{
  char* name;
  void* context;
  st_func_vtable vtable;
  int num_comp;
  bool homogeneous;
  bool constant;

  // Functions for computing derivatives (up to 4).
  st_func_t* derivs[4];
};

static void st_func_free(void* ctx, void* dummy)
{
  st_func_t* func = ctx;
  if ((func->vtable.dtor != NULL) && (func->context != NULL))
    func->vtable.dtor(func->context);
  func->context = NULL;
  polymec_free(func->name);
}

st_func_t* st_func_new(const char* name, void* context, st_func_vtable vtable,
                       st_func_homogeneity_t homogeneity,
                       st_func_constancy_t constancy,
                       int num_comp)
{
  ASSERT(context != NULL);
  ASSERT(vtable.eval != NULL);
  ASSERT(num_comp > 0);
  st_func_t* f = GC_MALLOC(sizeof(st_func_t));
  f->name = string_dup(name);
  f->context = context;
  f->vtable = vtable;
  f->homogeneous = (homogeneity == ST_HOMOGENEOUS);
  f->constant = (constancy == ST_CONSTANT);
  f->num_comp = num_comp;
  memset(f->derivs, 0, sizeof(st_func_t*)*4);
  GC_register_finalizer(f, st_func_free, f, NULL, NULL);
  return f;
}

st_func_t* st_func_from_func(const char* name, st_eval_func func, 
                             st_func_homogeneity_t homogeneity,
                             st_func_constancy_t constancy,
                             int num_comp)
{
  ASSERT(func != NULL);
  ASSERT(num_comp > 0);
  st_func_t* f = GC_MALLOC(sizeof(st_func_t));
  f->name = string_dup(name);
  f->context = NULL;
  f->vtable.eval = func;
  f->homogeneous = (homogeneity == ST_HOMOGENEOUS);
  f->constant = (constancy == ST_CONSTANT);
  f->num_comp = num_comp;
  memset(f->derivs, 0, sizeof(st_func_t*)*4);
  GC_register_finalizer(f, &st_func_free, f, NULL, NULL);
  return f;
}

static void eval_sp_func(void* context, point_t* x , real_t t, real_t* val)
{
  sp_func_t* sp_func = context;
  sp_func_eval(sp_func, x, val);
}

// This is used for converting sp_func derivatives to st_func derivatives. 
typedef struct
{
  sp_func_t* func;
  int deriv;
} st_func_sp_deriv_t;

static void st_func_sp_deriv_eval(void* context, point_t* x, real_t t, real_t* result)
{
  st_func_sp_deriv_t* F = context;
  sp_func_eval_deriv(F->func, F->deriv, x, result);
}

static void st_func_sp_deriv_free(void* context)
{
  st_func_sp_deriv_t* F = context;
  F->func = NULL;
  polymec_free(F);
}

static st_func_t* sp_func_deriv(sp_func_t* func, int d)
{
  ASSERT(d > 0);
  ASSERT(sp_func_has_deriv(func, d));
  st_func_sp_deriv_t* F = polymec_malloc(sizeof(st_func_sp_deriv_t));
  F->func = func;
  F->deriv = d;
  st_func_vtable vtable = {.eval = st_func_sp_deriv_eval, .dtor = st_func_sp_deriv_free};
  char name[1024];
  snprintf(name, 1024, "deriv(%s, %d)", sp_func_name(func), d);
  st_func_homogeneity_t homo = (sp_func_is_homogeneous(func)) ? ST_HOMOGENEOUS : ST_INHOMOGENEOUS;
  return st_func_new(name, F, vtable, homo, ST_CONSTANT, 3*sp_func_num_comp(func));
}

st_func_t* st_func_from_sp_func(sp_func_t* func)
{
  ASSERT(func != NULL);
  st_func_t* f = GC_MALLOC(sizeof(st_func_t));
  f->name = string_dup(sp_func_name(func));
  f->context = func;
  f->vtable.eval = eval_sp_func;
  f->homogeneous = sp_func_is_homogeneous(func);
  f->constant = true;
  f->num_comp = sp_func_num_comp(func);
  memset(f->derivs, 0, sizeof(st_func_t*)*4);
  for (int i = 1; i <= 4; ++i)
  {
    if (sp_func_has_deriv(func, i))
      st_func_register_deriv(f, i, sp_func_deriv(func, i));
  }
  GC_register_finalizer(f, &st_func_free, f, NULL, NULL);
  return f;
}

const char* st_func_name(st_func_t* func)
{
  return (const char*)func->name;
}

void st_func_rename(st_func_t* func, const char* new_name)
{
  polymec_free(func->name);
  func->name = string_dup(new_name);
}

bool st_func_is_homogeneous(st_func_t* func)
{
  return func->homogeneous;
}

bool st_func_is_constant(st_func_t* func)
{
  return func->constant;
}

int st_func_num_comp(st_func_t* func)
{
  return func->num_comp;
}

void* st_func_context(st_func_t* func)
{
  return func->context;
}

void st_func_eval(st_func_t* func, point_t* x, real_t t, real_t* result)
{
  func->vtable.eval(func->context, x, t, result);
}

void st_func_register_deriv(st_func_t* func, int n, st_func_t* nth_deriv)
{
  ASSERT(n > 0);
  ASSERT(n <= 4);
  ASSERT(nth_deriv != NULL);
  ASSERT(st_func_num_comp(nth_deriv) == (func->num_comp * (int)pow(3, n))); 
  func->derivs[n-1] = nth_deriv;
}

bool st_func_has_deriv(st_func_t* func, int n)
{
  ASSERT(n > 0);
  if (n > 4) return false; // FIXME
  return (func->derivs[n-1] != NULL);
}

// Evaluates the nth derivative of this function, placing the result in result.
void st_func_eval_deriv(st_func_t* func, int n, point_t* x, real_t t, real_t* result)
{
  st_func_eval(func->derivs[n-1], x, t, result);
}

typedef struct
{
  st_func_t* f; // Borrowed ref
  real_t t;
} st_frozen_ctx;

static void st_frozen_eval(void* ctx, point_t* x, real_t* result)
{
  st_frozen_ctx* c = (st_frozen_ctx*)ctx;
  st_func_eval(c->f, x, c->t, result); 
}

static void st_frozen_dtor(void* ctx)
{
  st_frozen_ctx* c = (st_frozen_ctx*)ctx;
  polymec_free(c);
}

sp_func_t* st_func_freeze(st_func_t* func, real_t t)
{
  sp_func_vtable vtable = {.eval = &st_frozen_eval, .dtor = &st_frozen_dtor };
  char name[1024];
  snprintf(name, 1024, "%s (frozen at %g)", st_func_name(func), t);
  st_frozen_ctx* c = polymec_malloc(sizeof(st_frozen_ctx));
  c->f = func;
  c->t = t;
  sp_func_homogeneity_t homog = (st_func_is_homogeneous(func)) ? SP_HOMOGENEOUS : SP_INHOMOGENEOUS;
  return sp_func_new(name, (void*)c, vtable, homog, st_func_num_comp(func));
}

// Multicomponent function stuff.

static int multicomp_st_func_magic_number = 23523462;

typedef struct 
{
  int magic_number;
  st_func_t** functions;
  int num_comp;
} multicomp_st_func_t;

static void multicomp_eval(void* context, point_t* x, real_t t, real_t* result)
{
  multicomp_st_func_t* mc = (multicomp_st_func_t*)context;
  for (int i = 0; i < mc->num_comp; ++i)
    st_func_eval(mc->functions[i], x, t, &result[i]);
}

static void multicomp_dtor(void* context)
{
  multicomp_st_func_t* mc = (multicomp_st_func_t*)context;
  polymec_free(mc->functions);
  polymec_free(mc);
}

st_func_t* multicomp_st_func_from_funcs(const char* name, 
                                        st_func_t** functions,
                                        int num_comp)
{
  ASSERT(num_comp >= 0);

  st_func_homogeneity_t homogeneity = ST_HOMOGENEOUS;
  st_func_constancy_t constancy = ST_CONSTANT;
  // The functions determine the constancy and homogeneity.
  multicomp_st_func_t* mc = polymec_malloc(sizeof(multicomp_st_func_t));
  mc->magic_number = multicomp_st_func_magic_number;
  mc->num_comp = num_comp;
  mc->functions = polymec_malloc(sizeof(st_func_t*)*num_comp);
  for (int i = 0; i < num_comp; ++i)
  {
    ASSERT(functions[i] != NULL);
    ASSERT(st_func_num_comp(functions[i]) == 1);
    mc->functions[i] = functions[i];
    if (!st_func_is_homogeneous(functions[i]))
      homogeneity = ST_INHOMOGENEOUS;
    if (!st_func_is_constant(functions[i]))
      constancy = ST_NONCONSTANT;
  }
  st_func_vtable vtable = {.eval = multicomp_eval, .dtor = multicomp_dtor};
  return st_func_new(name, (void*)mc, vtable, 
                     homogeneity, constancy, num_comp);
}

typedef struct 
{
  st_func_t* func;
  int num_comp, comp;
} extractedcomp_st_func_t;

static void extractedcomp_eval(void* context, point_t* x, real_t t, real_t* result)
{
  extractedcomp_st_func_t* ec = (extractedcomp_st_func_t*)context;
  real_t vals[ec->num_comp];
  st_func_eval(ec->func, x, t, vals);
  *result = vals[ec->comp];
}

static void extractedcomp_dtor(void* context)
{
  extractedcomp_st_func_t* ec = (extractedcomp_st_func_t*)context;
  polymec_free(ec);
}

st_func_t* st_func_from_component(st_func_t* multicomp_func,
                                  int component)
{
  ASSERT(component >= 0);
  ASSERT(component < st_func_num_comp(multicomp_func));
  
  // is this a proper multicomp_st_func_t?
  multicomp_st_func_t* mc = (multicomp_st_func_t*)multicomp_func;
  if (mc->magic_number == multicomp_st_func_magic_number)
  {
    // We can just use the original function directly!
    return mc->functions[component];
  }
  else
  {
    // We'll have to wrap this st_func.
    st_func_vtable vtable = {.eval = extractedcomp_eval, .dtor = extractedcomp_dtor};
    int name_len = strlen(st_func_name(multicomp_func)) + 10;
    char name[name_len];
    st_func_homogeneity_t homogeneity = st_func_is_homogeneous(multicomp_func) ? ST_HOMOGENEOUS : ST_INHOMOGENEOUS;
    st_func_constancy_t constancy = st_func_is_constant(multicomp_func) ? ST_CONSTANT : ST_NONCONSTANT;
    snprintf(name, name_len, "%s[%d]", st_func_name(multicomp_func), component);
    extractedcomp_st_func_t* ec = polymec_malloc(sizeof(extractedcomp_st_func_t));
    ec->func = multicomp_func;
    ec->num_comp = st_func_num_comp(multicomp_func);
    ec->comp = component;
    return st_func_new(name, (void*)ec, vtable, 
                       homogeneity, constancy, 1);
  }
}

// Constant space-time function implementation.

typedef struct
{
  int num_comp;
  real_t *comp;
} const_st_func_t;

static void constant_eval(void* ctx, point_t* x, real_t t, real_t* res)
{
  const_st_func_t* f = (const_st_func_t*)ctx;
  for (int i = 0; i < f->num_comp; ++i)
    res[i] = f->comp[i];
}

static void constant_dtor(void* ctx)
{
  const_st_func_t* f = (const_st_func_t*)ctx;
  polymec_free(f->comp);
  polymec_free(f);
}

static st_func_t* create_constant_st_func(real_t components[], int num_components)
{
  st_func_vtable vtable = {.eval = constant_eval, .dtor = constant_dtor};
  char name[1024];
  snprintf(name, 1024, "constant space-time function (");
  for (int i = 0; i < num_components; ++i)
  {
    char comp_str[22];
    if (i == (num_components-1))
      snprintf(comp_str, 20, "%g)", components[i]);
    else
      snprintf(comp_str, 20, "%g, ", components[i]);
    strncat(name, comp_str, 1024);
  }
  const_st_func_t* f = polymec_malloc(sizeof(const_st_func_t));
  f->num_comp = num_components;
  f->comp = polymec_malloc(num_components*sizeof(real_t));
  for (int i = 0; i < num_components; ++i)
    f->comp[i] = components[i];
  return st_func_new(name, (void*)f, vtable, ST_HOMOGENEOUS, ST_CONSTANT, num_components);
}

st_func_t* constant_st_func_new(real_t components[], int num_components)
{
  st_func_t* func = create_constant_st_func(components, num_components);

  // Just to be complete, we register the 3-component zero function as this 
  // function's gradient.
  real_t zeros[3*num_components];
  memset(zeros, 0, 3*num_components*sizeof(real_t));
  st_func_t* zero = create_constant_st_func(zeros, 3*num_components);
  st_func_register_deriv(func, 1, zero);

  return func;
}


