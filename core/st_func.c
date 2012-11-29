#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

struct st_func_t 
{
  char* name;
  void* context;
  st_vtable vtable;
  int num_comp;
  bool homogeneous;
  bool constant;
};

static void st_func_free(void* ctx, void* dummy)
{
  st_func_t* func = (st_func_t*)ctx;
  if (func->vtable.dtor)
    func->vtable.dtor(func->context);
  free(func->name);
  free(func);
}

st_func_t* st_func_new(const char* name, void* context, st_vtable vtable,
                       st_func_homogeneity_t homogeneity,
                       st_func_constancy_t constancy,
                       int num_comp)
{
  ASSERT(context != NULL);
  ASSERT(vtable.eval != NULL);
  ASSERT(num_comp > 0);
  st_func_t* f = GC_MALLOC(sizeof(st_func_t));
  f->name = strdup(name);
  f->context = context;
  f->vtable = vtable;
  f->homogeneous = (homogeneity == ST_HOMOGENEOUS);
  f->constant = (constancy == ST_CONSTANT);
  f->num_comp = num_comp;
  GC_register_finalizer(f, &st_func_free, f, NULL, NULL);
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
  f->name = strdup(name);
  f->context = NULL;
  f->vtable.eval = func;
  f->homogeneous = (homogeneity == ST_HOMOGENEOUS);
  f->constant = (constancy == ST_CONSTANT);
  f->num_comp = num_comp;
  GC_register_finalizer(f, &st_func_free, f, NULL, NULL);
  return f;
}

const char* st_func_name(st_func_t* func)
{
  return (const char*)func->name;
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

void st_func_eval(st_func_t* func, point_t* x, double t, double* result)
{
  func->vtable.eval(func->context, x, t, result);
}

typedef struct
{
  st_func_t* f; // Borrowed ref
  double t;
} st_frozen_ctx;

static void st_frozen_eval(void* ctx, point_t* x, double* result)
{
  st_frozen_ctx* c = (st_frozen_ctx*)ctx;
  st_func_eval(c->f, x, c->t, result); 
}

static void st_frozen_dtor(void* ctx)
{
  st_frozen_ctx* c = (st_frozen_ctx*)ctx;
  free(c);
}

sp_func_t* st_func_freeze(st_func_t* func, double t)
{
  sp_vtable vtable = {.eval = &st_frozen_eval, .dtor = &st_frozen_dtor };
  char name[1024];
  snprintf(name, 1024, "%s (frozen at %g)", st_func_name(func), t);
  st_frozen_ctx* c = malloc(sizeof(st_frozen_ctx));
  c->f = func;
  c->t = t;
  sp_func_homogeneity_t homog = (st_func_is_homogeneous(func)) ? SP_HOMOGENEOUS : SP_INHOMOGENEOUS;
  return sp_func_new(name, (void*)c, vtable, homog, st_func_num_comp(func));
}

// Multicomponent function stuff.

typedef struct 
{
  st_func_t** functions;
  int num_comp;
} multicomp_st_func_t;

static void multicomp_eval(void* context, point_t* x, double t, double* result)
{
  multicomp_st_func_t* mc = (multicomp_st_func_t*)context;
  for (int i = 0; i < mc->num_comp; ++i)
    st_func_eval(mc->functions[i], x, t, &result[i]);
}

static void multicomp_dtor(void* context)
{
  multicomp_st_func_t* mc = (multicomp_st_func_t*)context;
  free(mc->functions);
  free(mc);
}

st_func_t* multicomp_st_func_from_funcs(const char* name, 
                                        st_func_t** functions,
                                        int num_comp)
{
  st_func_homogeneity_t homogeneity = ST_HOMOGENEOUS;
  st_func_constancy_t constancy = ST_CONSTANT;
  // The functions determine the constancy and homogeneity.
  multicomp_st_func_t* mc = malloc(sizeof(multicomp_st_func_t));
  mc->num_comp = num_comp;
  mc->functions = malloc(sizeof(st_func_t*)*num_comp);
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
  st_vtable vtable = {.eval = multicomp_eval, .dtor = multicomp_dtor};
  return st_func_new(name, (void*)mc, vtable, 
                     homogeneity, constancy, num_comp);
}

#ifdef __cplusplus
}
#endif


