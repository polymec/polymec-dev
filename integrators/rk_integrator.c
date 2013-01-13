#include "integrators/rk_integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  void* context;
  rk_compute_deriv compute_deriv;
  integrator_dtor dtor;
  double* k;
} rk1_t;

static void rk1_dtor(void* context)
{
  rk1_t* rk1 = context;
  if ((rk1->context != NULL) && (rk1->dtor != NULL))
    rk1->dtor(rk1->context);
  free(rk1->k);
  free(rk1);
}

static void rk1_init(void* context, double t, double* solution, int N)
{
  rk1_t* rk1 = context;
  rk1->k = malloc(sizeof(double)*N);
}

static void rk1_step(void* context, double t1, double t2, double* solution, int N)
{
  ASSERT(t2 > t1);
  rk1_t* rk1 = context;
  rk1->compute_deriv(rk1->context, t1, solution, rk1->k);
  for (int i = 0; i < N; ++i)
    solution[N] += rk1->k[N];
}

static integrator_t* rk1_integrator_new(void* context, 
                                        rk_compute_deriv compute_deriv,
                                        integrator_dtor dtor)
{
  ASSERT(compute_deriv != NULL);
  rk1_t* rk1 = malloc(sizeof(rk1_t));
  rk1->context = context;
  rk1->compute_deriv = compute_deriv;
  rk1->dtor = dtor;
  integrator_vtable vtable = {.init = rk1_init, .step = rk1_step, .dtor = rk1_dtor};
  return integrator_new("RK1", rk1, vtable, 1, INTEGRATOR_EXPLICIT);
}

integrator_t* rk_integrator_new(int order,
                                void* context, 
                                rk_compute_deriv compute_deriv,
                                integrator_dtor dtor)
{
  ASSERT(order >= 1);
  ASSERT(order <= 4);
  switch (order)
  {
    case 1: return rk1_integrator_new(context, compute_deriv, dtor);
    default: return NULL;
  }
}

#ifdef __cplusplus
}
#endif

