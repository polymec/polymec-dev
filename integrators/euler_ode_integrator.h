// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_EULER_ODE_INTEGRATOR_H
#define POLYMEC_EULER_ODE_INTEGRATOR_H

#include "integrators/newton_solver.h"
#include "integrators/ode_integrator.h"

// This type of ODE integrator integrates a set of ordinary 
// differential equations using an implicit Euler method.

// Creates an Euler integrator that uses functional, or fixed-point, iteration. 
// The implicitness parameter theta determines the degree of implicitness, with 
// 0 signifying a forward Euler method, 0.5 a Crank-Nicolson method, and 
// 1 a backward Euler method. The integrator is constructed for a system of 
// ordinary differential equations (with given numbers of local and remote 
// values, for parallel capability), using a context pointer, a 
// right-hand side function, and a destructor. 
ode_integrator_t* functional_euler_ode_integrator_new(real_t theta,
                                                      MPI_Comm comm,
                                                      int num_local_values,
                                                      int num_remote_values,
                                                      void* context, 
                                                      int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                      void (*dtor)(void* context));

// Creates a backward Euler integrator that uses an inexact Newton-Krylov 
// solver with the requested timestep.
ode_integrator_t* newton_euler_ode_integrator_new(MPI_Comm comm,
                                                  int num_local_values,
                                                  int num_remote_values,
                                                  void* context,
                                                  int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                  void (*dtor)(void* context),
                                                  newton_pc_t* precond,
                                                  newton_krylov_t solver_type,
                                                  int max_krylov_dim);

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
void euler_ode_integrator_set_lp_convergence_norm(ode_integrator_t* integrator,
                                                  int p);

// Sets a custom norm to be used for measuring the change in the solution 
// over an iteration. The norm function takes an MPI communicator (over which 
// the norm should be summed/maxed/etc), vectors x and y to be differenced, 
// the length N of the vectors, and pointers to storage for the absolute and 
// relative norms.
void euler_ode_integrator_set_custom_convergence_norm(ode_integrator_t* integrator,
                                                      void (*compute_norms)(MPI_Comm comm, real_t* x, real_t* y, int N, real_t* abs_norm, real_t* rel_norm));

// This method provides access to the underlying Newton solver in the Newton Euler 
// ODE integrator.
newton_solver_t* newton_euler_ode_integrator_solver(ode_integrator_t* integrator);

// This observer type can be used to define objects that respond to actions
// taken by the euler_ode_integrator.
typedef struct euler_ode_observer_t euler_ode_observer_t;

// Creates and returns a newly-allocated observer that observes an 
// am_ode_integrator. Responses are given by the following arguments:
// rhs_computed - This function is called when the right hand side of the ODE
//                system is computed by the integrator.
// These functions are fed the given context object.
euler_ode_observer_t* euler_ode_observer_new(void* context,
                                             void (*rhs_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                             void (*dtor)(void* context));

// Destroys the given observer.
void euler_ode_observer_free(euler_ode_observer_t* observer);

// Adds the given observer to the given euler_ode_integrator. The observer 
// is consumed.
void euler_ode_integrator_add_observer(ode_integrator_t* integrator,
                                       euler_ode_observer_t* observer);

#endif

