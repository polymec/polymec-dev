// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_EULER_ODE_SOLVER_H
#define POLYMEC_EULER_ODE_SOLVER_H

#include "solvers/newton_solver.h"
#include "solvers/ode_solver.h"

// This type of ODE solver integrates a set of ordinary 
// differential equations using an implicit Euler method.

// Creates an Euler solver that uses functional, or fixed-point, iteration. 
// The implicitness parameter theta determines the degree of implicitness, with 
// 0 signifying a forward Euler method, 0.5 a Crank-Nicolson method, and 
// 1 a backward Euler method. The solver is constructed for a system of 
// ordinary differential equations (with given numbers of local and remote 
// values, for parallel capability), using a context pointer, a 
// right-hand side function, and a destructor. 
ode_solver_t* functional_euler_ode_solver_new(real_t theta,
                                              MPI_Comm comm,
                                              int num_local_values,
                                              int num_remote_values,
                                              void* context, 
                                              int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                              void (*dtor)(void* context));

// Creates a backward Euler solver that uses an Jacobian-Free Newton-Krylov 
// solver with the requested timestep.
ode_solver_t* jfnk_euler_ode_solver_new(MPI_Comm comm,
                                        int num_local_values,
                                        int num_remote_values,
                                        void* context,
                                        int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                        void (*dtor)(void* context),
                                        newton_pc_t* precond,
                                        jfnk_newton_t solver_type,
                                        int max_krylov_dim);

// Sets the maximum number of iterations in an integration step.
// The default is 100.
void euler_ode_solver_set_max_iterations(ode_solver_t* solver,
                                         int max_iters);

// Sets the relative and absolute tolerances that determine successful
// convergence of the iteration. The defaults are 1e-4 and 1.0, respectively.
void euler_ode_solver_set_tolerances(ode_solver_t* solver,
                                     real_t relative_tol, 
                                     real_t absolute_tol);

// Sets the Lp norm to be used for measuring the change in the solution over 
// an iteration. By default, the L-infinity norm is used. p may be 1 or 2, or 
// 0 (L-infinity norm).
void euler_ode_solver_set_lp_convergence_norm(ode_solver_t* solver,
                                              int p);

// Sets a custom norm to be used for measuring the change in the solution 
// over an iteration. The norm function takes an MPI communicator (over which 
// the norm should be summed/maxed/etc), vectors x and y to be differenced, 
// the length N of the vectors, and pointers to storage for the absolute and 
// relative norms.
void euler_ode_solver_set_custom_convergence_norm(ode_solver_t* solver,
                                                  void (*compute_norms)(MPI_Comm comm, real_t* x, real_t* y, int N, real_t* abs_norm, real_t* rel_norm));

// This method provides access to the underlying Newton solver in the Newton Euler 
// ODE solver.
newton_solver_t* newton_euler_ode_solver_solver(ode_solver_t* solver);

// This observer type can be used to define objects that respond to actions
// taken by the euler_ode_solver.
typedef struct euler_ode_observer_t euler_ode_observer_t;

// Creates and returns a newly-allocated observer that observes an 
// am_ode_solver. Responses are given by the following arguments:
// rhs_computed - This function is called when the right hand side of the ODE
//                system is computed by the solver.
// These functions are fed the given context object.
euler_ode_observer_t* euler_ode_observer_new(void* context,
                                             void (*rhs_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                             void (*dtor)(void* context));

// Destroys the given observer.
void euler_ode_observer_free(euler_ode_observer_t* observer);

// Adds the given observer to the given euler_ode_solver. The observer 
// is consumed.
void euler_ode_solver_add_observer(ode_solver_t* solver,
                                   euler_ode_observer_t* observer);

#endif

