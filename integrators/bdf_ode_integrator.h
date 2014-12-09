// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_BDF_ODE_INTEGRATOR_H
#define POLYMEC_BDF_ODE_INTEGRATOR_H

#include "integrators/ode_integrator.h"
#include "core/preconditioner.h"

// This type of ODE integrator integrates a stiff set of ordinary 
// differential equations using Backwards Difference Formulae (BDF). These 
// formulae can provide accurate time integration at an order from 1 to 5, and 
// L stability at orders 1 and 2. The BDF integrators all use Newton iteration.
// These integrators are implemented using CVODE from the Sundials suite 
// of nonlinear solvers.

// This function creates a BDF integrator that uses a Jacobian-Free
// Newton-Krylov method to solve the underlying linearized equations. This 
// method requires a preconditioner that is a coarse approximation of the 
// Jacobian matrix, but captures its essential behavior. The right-hand side 
// function and destructor are required, and a function to compute Jy, the 
// product of the Jacobian with the vector y, can be optionally provided. Its 
// signature is: 
// Jy(context, t, x, rhstx, y, temp, Jy), where context is the context pointer, 
// t is the time at which Jy is evaluated, y is the vector the Jacobian operator
// J is applied to, rhsty is the right hand side function evaluated at t and y,
// temp is a work vector the same size as y, and Jy is a vector that stores the 
// product Jy.
// If Jy is not given, a finite difference approximation of Jy will be used.
// Additionally, the type of Krylov solver (JFNK_BDF_GMRES, JFNK_BDF_BICGSTAB, 
// or JFNK_BDF_TFQMR) must be given, along with the maximum dimension of the 
// Krylov subspace. 
typedef enum 
{
  JFNK_BDF_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_BDF_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_BDF_TFQMR     // Transpose-Free QMR Krylov solver
} jfnk_bdf_krylov_t;

ode_integrator_t* jfnk_bdf_ode_integrator_new(int order, 
                                              MPI_Comm comm,
                                              int num_local_values, 
                                              int num_remote_values, 
                                              void* context, 
                                              int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                              int (*Jy)(void* context, real_t t, real_t* x, real_t* rhstx, real_t* y, real_t* temp, real_t* Jy),
                                              void (*dtor)(void* context),
                                              preconditioner_t* precond,
                                              jfnk_bdf_krylov_t solver_type,
                                              int max_krylov_dim);

// This returns the context pointer passed to the bdf_ode_integrator 
// constructor. In general, this will NOT return the same pointer as 
// ode_integrator_context!
void* bdf_ode_integrator_context(ode_integrator_t* integrator);

// This function evaluates error weights for use in the WRMS error norm.
typedef void (*bdf_ode_integrator_error_weight_func)(void* context, real_t* y, real_t* weights);

// Sets the relative and absolute tolerances for integrated quantities.
void bdf_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                       real_t relative_tol, real_t absolute_tol);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void bdf_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                  bdf_ode_integrator_error_weight_func compute_weights);                               

// Toggles stability limit detection, making the BDF method more robust at 
// a cost of 10-20% overhead.
void bdf_ode_integrator_set_stability_limit_detection(ode_integrator_t* integrator,
                                                      bool use_detection);

// Sets the maximum number of error test failures permitted in attempting 
// a single time step. By default, this value is 7.
void bdf_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                  int max_failures);

// Sets the maximum number of nonlinear solver iterations per time step.
// By default, this value is 3.
void bdf_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                     int max_iterations);

// Sets the safety factor (coefficient) used in the nonlinear convergence test.
// By default, this value is 0.1.
void bdf_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                        real_t coefficient);

// Evaluates the right-hand side of the system at the given time and with the 
// given solution X, placing the results in rhs.
void bdf_ode_integrator_eval_rhs(ode_integrator_t* integ, real_t t, real_t* X, real_t* rhs);

// Diagnostics for the time integrator.
typedef struct
{
  char* status_message; // borrowed pointer from integrator: do not free.
  long int num_steps;
  int order_of_last_step;
  real_t last_step_size;
  long int num_rhs_evaluations;
  long int num_linear_solve_setups;
  long int num_linear_solve_iterations;
  long int num_linear_solve_convergence_failures;
  long int num_error_test_failures;
  long int num_nonlinear_solve_iterations;
  long int num_nonlinear_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
} bdf_ode_integrator_diagnostics_t;

// Retrieve diagnostics for the time integrator.
void bdf_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                        bdf_ode_integrator_diagnostics_t* diagnostics);

// Writes time integrator diagnostics to the given file.
void bdf_ode_integrator_diagnostics_fprintf(bdf_ode_integrator_diagnostics_t* diagnostics, 
                                            FILE* stream);

#endif

