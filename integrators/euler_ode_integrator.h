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
// right-hand side function, and a destructor. 
ode_integrator_t* functional_euler_ode_integrator_new(real_t alpha,
                                                      MPI_Comm comm,
                                                      int num_local_values,
                                                      int num_remote_values,
                                                      void* context, 
                                                      int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                      void (*dtor)(void* context));

// Sets the maximum number of iterations in an integration step.
// The default is 100.
void euler_ode_integrator_set_max_iterations(ode_integrator_t* integrator,
                                             int max_iters);

// Sets the relative and absolute tolerances that determine successful
// convergence of the iteration. The defaults are 1e-4 and 1.0, respectively.
void euler_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                         real_t relative_tol, 
                                         real_t absolute_tol);

// Sets the Lp norm to be used for measuring the change in the solution over 
// an iteration. By default, the L-infinity norm is used. p may be 1 or 2, or 
// 0 (L-infinity norm).
void euler_ode_integrator_set_convergence_norm(ode_integrator_t* integrator,
                                               int p);

// By default, the Euler ODE integrator makes no attempt to adjust the 
// given timestep size. If this is called with a factor between 0 and 1, 
// the timestep will be adjusted to this factor times its old value when 
// an iteration fails.
void euler_ode_integrator_enable_dt_cuts(ode_integrator_t* integrator,
                                         real_t dt_cut_factor);
#endif

