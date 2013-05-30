#ifndef POLYMEC_RK_INTEGRATOR_H
#define POLYMEC_RK_INTEGRATOR_H

#include "integrators/integrator.h"

// A function signature for explicitly computing the derivative of a solution 
// at the given time. Arguments are:
// 1. A context object
// 2. The time t
// 3. The solution u 
// 4. Storage for the derivative of the solution, du/dt.
typedef void (*rk_compute_deriv)(void*, double, double*, double*);

// Creates an explicit Runge-Kutta integrator of the given order.
integrator_t* rk_integrator_new(int order,
                                void* context, 
                                rk_compute_deriv compute_deriv,
                                integrator_dtor dtor);

#endif

