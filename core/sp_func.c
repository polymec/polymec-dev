#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sp_func_t 
{
  char* name;
  void* context;
  sp_vtable vtable;
  int num_comp;
  bool homogeneous;
};

static void sp_func_free(void* ctx, void* dummy)
{
  UNUSED_ARG(dummy);
  sp_func_t* func = (sp_func_t*)ctx;
  if (func->vtable.dtor)
    free(func->context);
  free(func->name);
  free(func);
}

sp_func_t* sp_func_new(const char* name, void* context, sp_vtable vtable,
                       sp_func_homogeneity_t homogeneity,
                       int num_comp)
{
  ASSERT(context != NULL);
  ASSERT(vtable.eval != NULL);
  ASSERT(num_comp > 0);
  sp_func_t* f = GC_MALLOC(sizeof(sp_func_t));
  f->name = strdup(name);
  f->context = context;
  f->vtable = vtable;
  f->homogeneous = (homogeneity == SP_HOMOGENEOUS);
  f->num_comp = num_comp;
  GC_register_finalizer(f, &sp_func_free, f, NULL, NULL);
  return f;
}

sp_func_t* sp_func_from_func(const char* name, sp_eval_func func, 
                             sp_func_homogeneity_t homogeneity,
                             int num_comp)
{
  ASSERT(func != NULL);
  ASSERT(num_comp > 0);
  sp_func_t* f = GC_MALLOC(sizeof(sp_func_t));
  f->name = strdup(name);
  f->context = NULL;
  f->vtable.eval = func;
  f->homogeneous = (homogeneity == SP_HOMOGENEOUS);
  f->num_comp = num_comp;
  GC_register_finalizer(f, &sp_func_free, f, NULL, NULL);
  return f;
}

const char* sp_func_name(sp_func_t* func)
{
  return (const char*)func->name;
}

bool sp_func_is_homogeneous(sp_func_t* func)
{
  return func->homogeneous;
}

int sp_func_num_comp(sp_func_t* func)
{
  return func->num_comp;
}

void sp_func_eval(sp_func_t* func, point_t* x, double* result)
{
  func->vtable.eval(func->context, x, result);
}

#ifdef __cplusplus
}
#endif


