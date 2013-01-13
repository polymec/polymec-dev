#ifndef POLYMEC_FEULER_INTEGRATOR_H
#define POLYMEC_FEULER_INTEGRATOR_H

#include "core/integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

// A function signature for explicitly computing the derivative of a solution 
// at the given time. Arguments are:
// 1. A context object
// 2. The time t
// 3. The solution u 
// 4. Storage for the derivative of the solution, du/dt.
typedef void (*feuler_compute_deriv)(void*, double, double*, double*);

// Creates an integrator that uses the 1st-order explicit forward Euler method,
// using the given function to compute the derivative of the solution at t1.
integrator_t* feuler_integrator_new(void* context, 
                                    feuler_compute_deriv compute_deriv,
                                    integrator_dtor dtor);

#ifdef __cplusplus
}
#endif

#endif

