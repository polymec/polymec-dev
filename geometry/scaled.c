#include "geometry/scaled.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  sp_func_t* func;
  double scale_factor;
} sc_t;

static void sc_eval(void* ctx, point_t* x, double* result)
{
  sc_t* sc = ctx;
  // Scale x.
  double A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sp_func_eval(sc->func, &y, result);
}

static void sc_eval_gradient(void* ctx, point_t* x, double* result)
{
  sc_t* sc = ctx;
  ASSERT(sp_func_has_deriv(sc->func, 1));

  // Scale x.
  double A = sc->scale_factor;
  point_t y = {.x = x->x/A, .y = x->y/A, .z = x->z/A};
  sp_func_eval_deriv(sc->func, 1, &y, result);
}

sp_func_t* scaled_new(sp_func_t* func, double scale_factor)
{
  ASSERT(func != NULL);
  ASSERT(sp_func_num_comp(func) == 1);
  ASSERT(scale_factor > 0.0);

  sc_t* sc = malloc(sizeof(sc_t));
  sc->func = func;
  sc->scale_factor = scale_factor;
  char sc_str[1024];
  sprintf(sc_str, "scaled"); // FIXME: Not very helpful.
  sp_vtable vtable = {.eval = sc_eval};
  sp_func_t* sc_func = sp_func_new(sc_str, sc, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (sp_func_has_deriv(func, 1))
  {
    char sc_grad_str[1024];
    sprintf(sc_grad_str, "scaled gradient"); // FIXME: Yadda
    sp_vtable vtable_g = {.eval = sc_eval_gradient}; // Notice no dtor.
    sp_func_t* sc_grad = sp_func_new(sc_grad_str, sc, vtable_g, SP_INHOMOGENEOUS, 3);
    sp_func_register_deriv(sc_func, 1, sc_grad);
  }

  return sc_func;
}

#ifdef __cplusplus
}
#endif

