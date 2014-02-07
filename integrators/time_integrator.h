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

#ifndef POLYMEC_TIME_INTEGRATOR_H
#define POLYMEC_TIME_INTEGRATOR_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "integrators/preconditioner.h"

// This function evaluates the RHS of a system of ODEs.
typedef int (*time_integrator_rhs_func)(void* context, real_t t, real_t* x, real_t* x_dot);

// This function performs any parallel communication that is needed to 
// put values in place for the evaluation of the residual function.
typedef void (*time_integrator_communication_func)(void* context, real_t t, real_t* x);

// This function destroys the state (context).
typedef void (*time_integrator_dtor)(void* context);

// This function evaluates error weights for use in the WRMS error norm.
typedef void (*time_integrator_error_weight_func)(void* context, real_t* y, real_t* weights);

// This virtual table determines the behavior of the time integrator.
typedef struct
{
  // This function evaluates the right hand side of a coupled system of 
  // nonlinear ordinary different equations at time t with solution x, 
  // storing it in x_dot. It should return 0 on success, 1 for a 
  // recoverable error, -1 for a fatal error.
  time_integrator_rhs_func rhs;

  // Perform inter-process communication.
  time_integrator_communication_func communicate;

  // This (optional) function destroys the state (context) when the time integrator 
  // is destroyed.
  time_integrator_dtor dtor;

} time_integrator_vtable;

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. 
typedef struct time_integrator_t time_integrator_t;

// Creates an integrator that uses a GMRES Krylov method with a given 
// maximum subspace dimension of max_krylov_dim. No restarts are used.
// N is the dimension of the system.
time_integrator_t* gmres_time_integrator_new(const char* name, 
                                             void* context,
                                             MPI_Comm comm,
                                             int N,
                                             time_integrator_vtable vtable,
                                             int order,
                                             int max_krylov_dim);

// Creates an integrator that uses a stabilized bi-conjugate gradient Krylov 
// method with a given maximum subspace dimension of max_krylov_dim.
// N is the dimension of the system.
time_integrator_t* bicgs_time_integrator_new(const char* name, 
                                             void* context,
                                             MPI_Comm comm,
                                             int N,
                                             time_integrator_vtable vtable,
                                             int order,
                                             int max_krylov_dim);

// Creates an integrator that uses a transpose-free quasi-minimum residual 
// Krylov method with a given maximum subspace dimension of max_krylov_dim.
// N is the dimension of the system.
time_integrator_t* tfqmr_time_integrator_new(const char* name, 
                                             void* context,
                                             MPI_Comm comm,
                                             int N,
                                             time_integrator_vtable vtable,
                                             int order,
                                             int max_krylov_dim);

// Frees a time integrator.
void time_integrator_free(time_integrator_t* integrator);

// Returns the name of the integrator (internally stored).
char* time_integrator_name(time_integrator_t* integrator);

// Returns the context object for this integrator.
void* time_integrator_context(time_integrator_t* integrator);

// Returns the order of the integration method.
int time_integrator_order(time_integrator_t* integrator);

// Sets the preconditioner to use to help solve the equations.
void time_integrator_set_preconditioner(time_integrator_t* integrator,
                                        preconditioner_t* precond);

// Returns the current preconditioner matrix for the integrator, or 
// NULL if there is no such matrix.
preconditioner_matrix_t* time_integrator_preconditioner_matrix(time_integrator_t* integrator);

// Sets whether to use a stability limit detection algorithm to improve 
// robustness on particularly stiff problems (2-10% overhead, depending 
// on the problem).
void time_integrator_set_stability_limit_detection(time_integrator_t* integrator,
                                                   bool use_detection);

// Sets the relative and absolute tolerances for integrated quantities.
void time_integrator_set_tolerances(time_integrator_t* integrator,
                                    real_t relative_tol, real_t absolute_tol);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void time_integrator_set_error_weight_function(time_integrator_t* integrator,
                                               time_integrator_error_weight_func compute_weights);                               

// Evaluates the right-hand side of the system at the given time and with the 
// given solution X, placing the results in rhs.
void time_integrator_eval_rhs(time_integrator_t* integ, real_t t, real_t* X, real_t* rhs);

// Sets the maximum time step size for the next integration step.
void time_integrator_set_max_dt(time_integrator_t* integ, real_t max_dt);

// Sets the time past which the integrator will not step.
void time_integrator_set_stop_time(time_integrator_t* integ, real_t stop_time);

// Integrates the given solution X in place, taking a single step starting at 
// time *t and storing the new time in *t as well. Returns true if the step 
// succeeded, false if it failed for some reason. If a step fails, both t 
// and X remain unchanged.
bool time_integrator_step(time_integrator_t* integrator, real_t* t, real_t* X);

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
} time_integrator_diagnostics_t;

// Retrieve diagnostics for the time integrator.
void time_integrator_get_diagnostics(time_integrator_t* integrator, 
                                     time_integrator_diagnostics_t* diagnostics);

// Writes time integrator diagnostics to the given file.
void time_integrator_diagnostics_fprintf(time_integrator_diagnostics_t* diagnostics, 
                                         FILE* stream);

#endif

