#include "integrators/feuler_integrator.h"
#include "integrators/rk_integrator.h"

integrator_t* feuler_integrator_new(void* context, 
                                    feuler_compute_deriv compute_deriv,
                                    integrator_dtor dtor)
{
  // The Forward Euler method is implemented by the 1st-order RK method.
  return rk_integrator_new(1, context, compute_deriv, dtor);
}

