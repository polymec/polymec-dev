#include <stdlib.h>
#include <string.h>
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

st_func_t* st_func_new(const char* name, void* context, st_vtable vtable,
                       st_func_homogeneity_t homogeneity,
                       st_func_constancy_t constancy,
                       int num_comp)
{
  ASSERT(context != NULL);
  ASSERT(vtable.eval != NULL);
  ASSERT(num_comp > 0);
  st_func_t* f = malloc(sizeof(st_func_t));
  f->name = strdup(name);
  f->context = context;
  f->vtable = vtable;
  f->homogeneous = (homogeneity == ST_HOMOGENEOUS);
  f->constant = (constancy == ST_CONSTANT);
  f->num_comp = num_comp;
  return f;
}

st_func_t* st_func_from_func(const char* name, st_eval_func func, 
                             st_func_homogeneity_t homogeneity,
                             st_func_constancy_t constancy,
                             int num_comp)
{
  ASSERT(func != NULL);
  ASSERT(num_comp > 0);
  st_func_t* f = malloc(sizeof(st_func_t));
  f->name = strdup(name);
  f->context = NULL;
  f->vtable.eval = func;
  f->homogeneous = (homogeneity == ST_HOMOGENEOUS);
  f->constant = (constancy == ST_CONSTANT);
  f->num_comp = num_comp;
  return f;
}

void st_func_free(st_func_t* func)
{
  if (func->vtable.dtor)
    free(func->context);
  free(func->name);
  free(func);
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

#ifdef __cplusplus
}
#endif


