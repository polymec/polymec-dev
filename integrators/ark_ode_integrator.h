// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ARK_ODE_INTEGRATOR_H
#define POLYMEC_ARK_ODE_INTEGRATOR_H

#include "integrators/ode_integrator.h"
#include "integrators/newton_pc.h"

// This type of ODE integrator integrates a multi-scale set of ordinary 
// differential equations using semi-implicit Additive Runge Kutte (ARK) 
// methods. The equations integrated take the form
//
// M * dx/dt = fe(t, x) + fi(t, x)
//
// where M is a "mass matrix" term, fe is a slowly-varying function that 
// can be integrated explicitly, and fi is a quickly-varying function that 
// should be integrated implicitly.

typedef enum 
{
  JFNK_BDF_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_BDF_FGMRES,   // Flexible GMRES Krylov solver
  JFNK_BDF_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_BDF_TFQMR,    // Transpose-Free QMR Krylov solver
  JFNK_BDF_PCG       // Preconditioned Conjugate Gradient Krylov solver
} jfnk_ark_krylov_t;

ode_integrator_t* explicit_ark_ode_integrator_new(int order, 
                                                  MPI_Comm comm,
                                                  int num_local_values, 
                                                  int num_remote_values, 
                                                  void* context, 
                                                  int (*fe_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                                  void (*dtor)(void* context));

ode_integrator_t* functional_ark_ode_integrator_new(int order, 
                                                    MPI_Comm comm,
                                                    int num_local_values, 
                                                    int num_remote_values, 
                                                    void* context, 
                                                    int (*fe_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                                    int (*fi_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                                    void (*dtor)(void* context));

ode_integrator_t* jfnk_ark_ode_integrator_new(int order, 
                                              MPI_Comm comm,
                                              int num_local_values, 
                                              int num_remote_values, 
                                              void* context, 
                                              int (*fe_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                              int (*fi_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                              int (*Jy_func)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* temp, real_t* Jy),
                                              void (*dtor)(void* context),
                                              newton_pc_t* precond,
                                              jfnk_ark_krylov_t solver_type,
                                              int max_krylov_dim);

// This returns the context pointer passed to the ark_ode_integrator 
// constructor. In general, this will NOT return the same pointer as 
// ode_integrator_context!
void* ark_ode_integrator_context(ode_integrator_t* integrator);

// Sets the relative and absolute tolerances for integrated quantities.
void ark_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                       real_t relative_tol, real_t absolute_tol);

// Sets the error weights for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances. Weights are copied into the integrator.
void ark_ode_integrator_set_error_weights(ode_integrator_t* integrator, real_t* weights);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void ark_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                  void (*compute_weights)(void* context, real_t* y, real_t* weights));

// Sets the maximum number of error test failures permitted in attempting 
// a single time step. By default, this value is 7.
void ark_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                  int max_failures);

// Sets the maximum number of nonlinear solver iterations per time step.
// By default, this value is 3.
void ark_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                     int max_iterations);

// Sets the safety factor (coefficient) used in the nonlinear convergence test.
// By default, this value is 0.1.
void ark_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                        real_t coefficient);

// Evaluates the right-hand side of the system at the given time and with the 
// given solution X, placing the results in rhs.
void ark_ode_integrator_eval_rhs(ode_integrator_t* integ, real_t t, real_t* X, real_t* rhs);

// Returns an internal pointer to the preconditioner passed to this 
// integrator during construction time.
newton_pc_t* ark_ode_integrator_preconditioner(ode_integrator_t* integrator);

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
} ark_ode_integrator_diagnostics_t;

// Retrieve diagnostics for the time integrator.
void ark_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                        ark_ode_integrator_diagnostics_t* diagnostics);

// Writes time integrator diagnostics to the given file.
void ark_ode_integrator_diagnostics_fprintf(ark_ode_integrator_diagnostics_t* diagnostics, 
                                            FILE* stream);

// This observer type can be used to define objects that respond to actions
// taken by the ark_ode_integrator.
typedef struct ark_ode_observer_t ark_ode_observer_t;

// Creates and returns a newly-allocated observer that observes an 
// am_ode_integrator. Responses are given by the following arguments:
// rhs_computed - This function is called when the right hand side of the ODE
//                system is computed by the integrator.
// Jy_computed - This function is called when the Jacobian-vector product J*y 
//               is computed by the integrator.
// Both of these functions are fed the given context object.
ark_ode_observer_t* ark_ode_observer_new(void* context,
                                         void (*fe_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                         void (*fi_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                         void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy),
                                         void (*dtor)(void* context));

// Adds the given observer to the given am_ode_integrator. The observer 
// is consumed by the integrator, so no destructor is needed.
void ark_ode_integrator_add_observer(ode_integrator_t* integrator,
                                     ark_ode_observer_t* observer);

#endif

