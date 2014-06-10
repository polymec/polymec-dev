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

#ifndef POLYMEC_NONLINEAR_INTEGRATOR_H
#define POLYMEC_NONLINEAR_INTEGRATOR_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "integrators/preconditioner.h"

// The different global strategies for the Newton iteration.
typedef enum
{
  LINE_SEARCH,
  NO_GLOBAL_STRATEGY
} nonlinear_integrator_strategy_t;

// This is a pointer to a residual function used with a nonlinear integrator.
// This function evaluates the residual function for the nonlinear system
// of equations at the time t using the solution vector x and placing the 
// result in F. It should return 0 on success, 1 for a recoverable error, 
// -1 for a fatal error.
typedef int (*nonlinear_integrator_residual_func)(void* context, real_t t, real_t* x, real_t* F);

// Scaling function.
typedef void (*nonlinear_integrator_scaling_func)(void* context, real_t* scale_factor);

// Constraints function.
typedef void (*nonlinear_integrator_constraints_func)(void* context, real_t* constraints);

// Context destructor.
typedef void (*nonlinear_integrator_dtor)(void* context);

// This virtual table determines the behavior of the nonlinear integrator.
typedef struct
{
  nonlinear_integrator_residual_func eval;

  // This (optional) function sets the "x-scaling vector," which contains the diagonal 
  // components of a matrix Dx such that the components of Dx * x all have 
  // roughly the same magnitude as F(x) approaches 0.
  // - F_scale: the diagonal components of a matrix Df such that the components 
  //            of Df * F(x) all have roughly the same magnitude as F(x) approaches 0.
  nonlinear_integrator_scaling_func set_x_scale;

  // This (optional) function sets the "F-scaling vector," which contains the diagonal 
  // components of a matrix Df such that the components of Df * F(x) all have 
  // roughly the same magnitude as F(x) approaches 0.
  nonlinear_integrator_scaling_func set_F_scale;

  // This (optional) function sets the contraints vector, which places algebraic 
  // constraints on the components of the solution vector x. If constraints[i] is:
  // 0.0  - no constraint is placed on x[i].
  // 1.0  - x[i] must be non-negative.
  // -1.0 - x[i] must be non-positive.
  // 2.0  - x[i] must be positive.
  nonlinear_integrator_constraints_func set_constraints;

  // This (optional) function destroys the state (context) when the nonlinear integrator 
  // is destroyed.
  nonlinear_integrator_dtor dtor;

} nonlinear_integrator_vtable;

// This class represents a collection of algorithms for integrating partial 
// differential equations that are discretized into a (sparse) system of 
// nonlinear equations. The integration is performed using Matrix-free 
// Newton-Krylov methods provided by KINSol.
typedef struct nonlinear_integrator_t nonlinear_integrator_t;

// The following functions each create a nonlinear integrator with the given 
// name, context, and virtual table for integrated a discretized system of 
// partial differential equations represented by a sparse nonlinear system of 
// equations. The virtual table defines accessor methods for the residual 
// function and the adjacency graph, and the globalization strategy can be 
// NONE or LINE_SEARCH.

// Creates an integrator that uses a GMRES Krylov method with a given 
// maximum subspace dimension of max_krylov_dim, and a maximum number of 
// restarts given by max_restarts. N is the dimension of the system.
nonlinear_integrator_t* gmres_nonlinear_integrator_new(const char* name,
                                                       void* context,
                                                       MPI_Comm comm,
                                                       int N,
                                                       nonlinear_integrator_vtable vtable,
                                                       nonlinear_integrator_strategy_t global_strategy,
                                                       int max_krylov_dim,
                                                       int max_restarts);

// Creates an integrator that uses a stabilized bi-conjugate gradient Krylov 
// method with a given maximum subspace dimension of max_krylov_dim.
// N is the dimension of the system.
nonlinear_integrator_t* bicgstab_nonlinear_integrator_new(const char* name,
                                                          void* context,
                                                          MPI_Comm comm,
                                                          int N,
                                                          nonlinear_integrator_vtable vtable,
                                                          nonlinear_integrator_strategy_t global_strategy,
                                                          int max_krylov_dim);

// Creates an integrator that uses a transpose-free quasi-minimum residual 
// Krylov method with a given maximum subspace dimension of max_krylov_dim.
// N is the dimension of the system.
nonlinear_integrator_t* tfqmr_nonlinear_integrator_new(const char* name,
                                                       void* context,
                                                       MPI_Comm comm,
                                                       int N,
                                                       nonlinear_integrator_vtable vtable,
                                                       nonlinear_integrator_strategy_t global_strategy,
                                                       int max_krylov_dim);

// Frees a solver.
void nonlinear_integrator_free(nonlinear_integrator_t* integrator);

// Returns an internal string storing the name of the solver.
char* nonlinear_integrator_name(nonlinear_integrator_t* integrator);

// Returns the context pointer for the solver.
void* nonlinear_integrator_context(nonlinear_integrator_t* integrator);

// Returns the number of (local) equations in the nonlinear system.
int nonlinear_integrator_num_equations(nonlinear_integrator_t* integrator);

// Sets the tolerances for the function norm (norm_tolerance) and the Newton
// step (step_tolerance) for the nonlinear integrator.
void nonlinear_integrator_set_tolerances(nonlinear_integrator_t* integrator, real_t norm_tolerance, real_t step_tolerance);

// Sets the maximum number of Newton iterations for the integrator.
void nonlinear_integrator_set_max_iterations(nonlinear_integrator_t* integrator, int max_iterations);

// Sets the preconditioner to use to help solve the nonlinear system.
void nonlinear_integrator_set_preconditioner(nonlinear_integrator_t* integrator,
                                             preconditioner_t* precond);

// Returns the current preconditioner matrix for the integrator, or 
// NULL if there is no such matrix.
preconditioner_matrix_t* nonlinear_integrator_preconditioner_matrix(nonlinear_integrator_t* integrator);

// Sets the null space associated with the nonlinear operator F(X).
// The null space consists of the set of vectors {X: X != 0}: dF/dX * X = 0.
// If the null space contains the spatially homogeneous functions, the 
// homogeneous_functions flag should be set to true, and no such homogeneous 
// vector should be included in the null space. null_dim is the dimension of 
// the null space excluding the homogeneous functions.
void nonlinear_integrator_set_null_space(nonlinear_integrator_t* integrator,
                                         bool homogeneous_functions,
                                         real_t** null_space_vectors,
                                         int null_dim);

// Evaluates the residual vector, storing it in F.
void nonlinear_integrator_eval_residual(nonlinear_integrator_t* integrator, real_t t, real_t* X, real_t* F);

// Integrates the nonlinear system of equations F(X, t) = 0 in place, 
// using X as the initial guess. Returns true if the solution was obtained, 
// false if not. The number of nonlinear iterations will be stored in 
// num_iterations upon success.
bool nonlinear_integrator_solve(nonlinear_integrator_t* integrator, real_t t, real_t* X, int* num_iterations);

// Diagnostics for the nonlinear integrator.
typedef struct
{
  char* status_message; // borrowed pointer from integrator: do not free.
  long int num_function_evaluations;
  long int num_beta_condition_failures;
  long int num_backtrack_operations;
  long int num_nonlinear_iterations;
  real_t scaled_function_norm;
  real_t scaled_newton_step_length;
  long int num_linear_solve_iterations;
  long int num_linear_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
  long int num_jacobian_vector_product_evaluations;
  long int num_difference_quotient_function_evaluations;
} nonlinear_integrator_diagnostics_t;

// Retrieve diagnostics for the nonlinear integrator.
void nonlinear_integrator_get_diagnostics(nonlinear_integrator_t* integrator, 
                                          nonlinear_integrator_diagnostics_t* diagnostics);

// Writes nonlinear integrator diagnostics to the given file.
void nonlinear_integrator_diagnostics_fprintf(nonlinear_integrator_diagnostics_t* diagnostics, 
                                              FILE* stream);

#endif

