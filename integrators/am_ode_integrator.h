// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_AM_ODE_INTEGRATOR_H
#define POLYMEC_AM_ODE_INTEGRATOR_H

#include "integrators/ode_integrator.h"
#include "integrators/newton_pc.h"

// This type of ODE integrator integrates a non-stiff set of ordinary 
// differential equations using the Adams-Moulton formulae. They can provide 
// accurate time integration at an order from 1 to 12. There are two kinds of 
// Adams-Moulton ODE integrators: those based on functional iteration, which 
// requires only function evaluations, and those based on Newton iteration, 
// which requires the representation of a Jacobian matrix. Both of these 
// types of integrators are implemented using CVODE from the Sundials suite 
// of nonlinear solvers.

// The functional iteration variant of the Adams-Moulton integrator is particularly 
// simple and can be constructed on the communicator comm for a number of 
// non-stiff ordinary differential equations (with given numbers of local and 
// remote values, for parallel capability), using a context pointer, a 
// right-hand side function, and a destructor.
ode_integrator_t* functional_am_ode_integrator_new(int order, 
                                                   MPI_Comm comm, 
                                                   int num_local_values,
                                                   int num_remote_values,
                                                   void* context, 
                                                   int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                   void (*dtor)(void* context));

// This function creates an Adams-Moulton integrator that uses a Jacobian-Free
// Newton-Krylov method to solve the underlying linearized equations. This 
// method requires a preconditioner that is a coarse approximation of the 
// Jacobian matrix, but captures its essential behavior. The right-hand side 
// function and destructor are required, and a function to compute Jy, the 
// product of the Jacobian with the vector y, can be optionally provided. Its 
// signature is: 
// Jy_func(context, t, x, rhs, y, temp, Jy), where context is the context pointer, 
// t is the time at which Jy is evaluated, y is the vector the Jacobian operator
// J is applied to, rhs is the right hand side function evaluated at t and y,
// temp is a work vector the same size as y, and Jy is a vector that stores the 
// product Jy.
// If Jy_func is not given, a finite difference approximation of Jy will be used.
// Additionally, the type of Krylov solver (JFNK_AM_GMRES, JFNK_AM_BICGSTAB, 
// or JFNK_AM_TFQMR) must be given, along with the maximum dimension of the 
// Krylov subspace. 
typedef enum 
{
  JFNK_AM_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_AM_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_AM_TFQMR     // Transpose-Free QMR Krylov solver
} jfnk_am_krylov_t;

ode_integrator_t* jfnk_am_ode_integrator_new(int order, 
                                             MPI_Comm comm,
                                             int num_local_equations, 
                                             int num_remote_values, 
                                             void* context, 
                                             int (*rhs_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                             int (*Jy_func)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* temp, real_t* Jy),
                                             void (*dtor)(void* context),
                                             newton_pc_t* precond,
                                             jfnk_am_krylov_t solver_type,
                                             int max_krylov_dim);

// This returns the context pointer passed to the am_ode_integrator 
// constructor. In general, this will NOT return the same pointer as 
// ode_integrator_context!
void* am_ode_integrator_context(ode_integrator_t* integrator);

// This function evaluates error weights for use in the WRMS error norm.
typedef void (*am_ode_integrator_error_weight_func)(void* context, real_t* y, real_t* weights);

// Sets the relative and absolute tolerances for integrated quantities.
void am_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                      real_t relative_tol, real_t absolute_tol);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void am_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                 am_ode_integrator_error_weight_func compute_weights);                               

// Sets the maximum number of error test failures permitted in attempting 
// a single time step. By default, this value is 7.
void am_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                 int max_failures);

// Sets the maximum number of nonlinear solver iterations per time step.
// By default, this value is 3.
void am_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                    int max_iterations);

// Sets the safety factor (coefficient) used in the nonlinear convergence test.
// By default, this value is 0.1.
void am_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                       real_t coefficient);

// Evaluates the right-hand side of the system at the given time and with the 
// given solution X, placing the results in rhs.
void am_ode_integrator_eval_rhs(ode_integrator_t* integ, real_t t, real_t* X, real_t* rhs);

// Returns an internal pointer to the preconditioner passed to this 
// integrator during construction time. (This returns NULL for the functional 
// integrator, which doesn't need a preconditioner.)
newton_pc_t* am_ode_integrator_preconditioner(ode_integrator_t* integrator);

// Diagnostics for the time integrator.
typedef struct
{
  char* status_message; // borrowed pointer from integrator: do not free.
  long int num_steps;
  int order_of_last_step, order_of_next_step;
  real_t last_step_size, next_step_size;
  long int num_rhs_evaluations;
  long int num_linear_solve_setups;
  long int num_linear_solve_iterations;
  long int num_linear_solve_convergence_failures;
  long int num_error_test_failures;
  long int num_nonlinear_solve_iterations;
  long int num_nonlinear_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
} am_ode_integrator_diagnostics_t;

// Retrieve diagnostics for the time integrator.
void am_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                       am_ode_integrator_diagnostics_t* diagnostics);

// Writes time integrator diagnostics to the given file.
void am_ode_integrator_diagnostics_fprintf(am_ode_integrator_diagnostics_t* diagnostics, 
                                           FILE* stream);

// This observer type can be used to define objects that respond to actions
// taken by the am_ode_integrator.
typedef struct am_ode_observer_t am_ode_observer_t;

// Creates and returns a newly-allocated observer that observes an 
// am_ode_integrator. Responses are given by the following arguments:
// rhs_computed - This function is called when the right hand side of the ODE
//                system is computed by the integrator.
// Jy_computed - This function is called when the Jacobian-vector product J*y 
//               is computed by the integrator.
// Both of these functions are fed the given context object.
am_ode_observer_t* am_ode_observer_new(void* context,
                                       void (*rhs_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                       void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy),
                                       void (*dtor)(void* context));

// Destroys the given observer.
void am_ode_observer_free(am_ode_observer_t* observer);

// Adds the given observer to the given am_ode_integrator. The observer 
// is consumed.
void am_ode_integrator_add_observer(ode_integrator_t* integrator,
                                    am_ode_observer_t* observer);

#endif

