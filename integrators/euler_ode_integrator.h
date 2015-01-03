// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_EULER_ODE_INTEGRATOR_H
#define POLYMEC_EULER_ODE_INTEGRATOR_H

#include "integrators/ode_integrator.h"

// This type of ODE integrator integrates a non-stiff set of ordinary 
// differential equations using an implicit Euler method.

// Creates an Euler integrator that uses functional, or fixed-point, iteration. 
// The implicitness parameter alpha determines the degree of implicitness, with 
// 0 signifying a forward Euler method, 0.5 a Crank-Nicolson method, and 
// 1 a backward Euler method. The integrator is constructed for a system of 
// ordinary differential equations (with given numbers of local and remote 
// values, for parallel capability), using a context pointer, a 
// right-hand side function, and a destructor. This integrator does not manage
// parallelism--it only leaves room for remote values that you specify, and 
// these values are assumed to be at the end of the solution vector.
ode_integrator_t* functional_euler_ode_integrator_new(real_t alpha, 
                                                      int num_local_values,
                                                      int num_remote_values,
                                                      void* context, 
                                                      int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                      void (*dtor)(void* context));

// Sets the maximum number of iterations in an integration step.
void euler_ode_integrator_set_max_iterations(ode_integrator_t* integrator,
                                             int max_iters);

// Sets the relative and absolute tolerances that determine successful
// convergence of the iteration.
void euler_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                         real_t relative_tol, 
                                         real_t absolute_tol);

#endif

