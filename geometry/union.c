#include "geometry/union.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  sp_func_t** funcs;
  int num_funcs;
} un_t;

// Destructor.
static void un_free(void* ctx)
{
  un_t* un = ctx;
  free(un->funcs);
  free(un);
}

static void un_eval(void* ctx, point_t* x, double* result)
{
  un_t* un = ctx;
  double minval = FLT_MAX;
  for (int i = 0; i < un->num_funcs; ++i)
  {
    double ival;
    sp_func_eval(un->funcs[i], x, &ival);
    minval = MIN(minval, ival);
  }
  result[0] = minval;
}

static void un_eval_gradient(void* ctx, point_t* x, double* result)
{
  un_t* un = ctx;
  double minval = FLT_MAX;
  int index = -1;
  for (int i = 0; i < un->num_funcs; ++i)
  {
    double ival;
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

  un_t* un = malloc(sizeof(un_t));
  un->funcs = malloc(sizeof(sp_func_t*) * num_surfaces);
  un->num_funcs = num_surfaces;
  bool has_grad = true;
  for (int i = 0; i < num_surfaces; ++i)
  {
    ASSERT(sp_func_num_comp(surfaces[i]) == 1);
    un->funcs[i] = surfaces[i];
    if (!sp_func_has_deriv(surfaces[i], 1))
      has_grad = false;
  }

  char un_str[1024];
  sprintf(un_str, "union"); // FIXME: Not very helpful.
  sp_vtable vtable = {.eval = un_eval, .dtor = un_free};
  sp_func_t* union_func = sp_func_new(un_str, un, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (has_grad)
  {
    char un_grad_str[1024];
    sprintf(un_grad_str, "union gradient"); // FIXME: Yadda
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

#ifdef __cplusplus
}
#endif

