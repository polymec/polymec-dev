#include "integrators/feuler_integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  void* context;
  feuler_compute_deriv compute_deriv;
  integrator_dtor dtor;
  double* deriv;
} feuler_t;

static void feuler_dtor(void* context)
{
  feuler_t* feuler = context;
  if ((feuler->context != NULL) && (feuler->dtor != NULL))
    feuler->dtor(feuler->context);
  free(feuler->deriv);
  free(feuler);
}

static void feuler_init(void* context, double t, double* solution, int N)
{
  feuler_t* feuler = context;
  feuler->deriv = malloc(sizeof(double)*N);
}

static void feuler_step(void* context, double t1, double t2, double* solution, int N)
{
  ASSERT(t2 > t1);
  feuler_t* feuler = context;
  feuler->compute_deriv(feuler->context, t1, solution, feuler->deriv);
  for (int i = 0; i < N; ++i)
    solution[N] += feuler->deriv[N];
}

integrator_t* feuler_integrator_new(void* context, 
                                    feuler_compute_deriv compute_deriv,
                                    integrator_dtor dtor)
{
  ASSERT(create_vector != NULL);
  ASSERT(axpy != NULL);
  ASSERT(compute_deriv != NULL);
  ASSERT(free_vector != NULL);
  feuler_t* feuler = malloc(sizeof(feuler_t));
  feuler->context = context;
  feuler->compute_deriv = compute_deriv;
  feuler->dtor = dtor;
  integrator_vtable vtable = {.init = feuler_init, .step = feuler_step, .dtor = feuler_dtor};
  return integrator_new("Forward Euler", feuler, vtable, 1, INTEGRATOR_EXPLICIT);
}

#ifdef __cplusplus
}
#endif

