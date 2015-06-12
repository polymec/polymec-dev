// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NEWTON_SOLVER_H
#define POLYMEC_NEWTON_SOLVER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "integrators/newton_pc.h"

// This is a pointer to a residual function used with a nonlinear solver.
// This function evaluates the residual function for the nonlinear system
// of equations at the time t using the solution vector x and placing the 
// result in F. It should return 0 on success, 1 for a recoverable error, 
// -1 for a fatal error.
typedef int (*newton_solver_residual_func)(void* context, real_t t, real_t* x, real_t* F);

// This is used to define the Jacobian-vector project. Should return 0 on 
// success, a positive value for a recoverable failure, and a negative number 
// for an unrecoverable failure.
typedef int (*newton_solver_Jv_func)(void* context, bool new_x, real_t t, real_t* x, real_t* v, real_t* Jv);

// Initial guess function.
typedef void (*newton_solver_initial_guess_func)(void* context, real_t t, real_t* x);

// Scaling function.
typedef void (*newton_solver_scaling_func)(void* context, real_t* scale_factor);

// Constraints function.
typedef void (*newton_solver_constraints_func)(void* context, real_t* constraints);

// Context destructor.
typedef void (*newton_solver_dtor)(void* context);

// This virtual table determines the behavior of the nonlinear solver.
typedef struct
{
  // This function evaluates the residual of the nonlinear system.
  newton_solver_residual_func eval;

  // This (optional) function sets the "x-scaling vector," which contains the diagonal 
  // components of a matrix Dx such that the components of Dx * x all have 
  // roughly the same magnitude as F(x) approaches 0.
  // - F_scale: the diagonal components of a matrix Df such that the components 
  //            of Df * F(x) all have roughly the same magnitude as F(x) approaches 0.
  newton_solver_scaling_func set_x_scale;

  // This (optional) function sets the "F-scaling vector," which contains the diagonal 
  // components of a matrix Df such that the components of Df * F(x) all have 
  // roughly the same magnitude as F(x) approaches 0.
  newton_solver_scaling_func set_F_scale;

  // This (optional) function sets the contraints vector, which places algebraic 
  // constraints on the components of the solution vector x. If constraints[i] is:
  // 0.0  - no constraint is placed on x[i].
  // 1.0  - x[i] must be non-negative.
  // -1.0 - x[i] must be non-positive.
  // 2.0  - x[i] must be positive.
  newton_solver_constraints_func set_constraints;

  // This (optional) function destroys the state (context) when the nonlinear solver 
  // is destroyed.
  newton_solver_dtor dtor;

} newton_solver_vtable;

// This class represents a collection of algorithms for integrating partial 
// differential equations that are discretized into a (sparse) system of 
// nonlinear equations. The integration is performed using Matrix-free 
// Newton-Krylov methods provided by KINSol with right preconditioning.
typedef struct newton_solver_t newton_solver_t;

typedef enum 
{
  NEWTON_GMRES,    // Generalized minimum residual Krylov solver
  NEWTON_FGMRES,   // Flexible GMRES Krylov solver
  NEWTON_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  NEWTON_TFQMR     // Transpose-Free QMR Krylov solver
} newton_krylov_t;

// Creates a Newton solver that uses a Newton-Krylov method with a given 
// maximum subspace dimension of max_krylov_dim. If the solver_type is NEWTON_GMRES
// or NEWTON_FGMRES, the maximum number of restarts is given by max_restarts--
// otherwise that parameter is ignored. 
newton_solver_t* newton_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   newton_solver_vtable vtable,
                                   newton_pc_t* precond,
                                   newton_krylov_t solver_type,
                                   int max_krylov_dim,
                                   int max_restarts);

// Frees a solver.
void newton_solver_free(newton_solver_t* solver);

// Returns the context pointer for the solver.
void* newton_solver_context(newton_solver_t* solver);

// Returns the number of (local) equations in the nonlinear system.
int newton_solver_num_equations(newton_solver_t* solver);

// Sets a function to use to evaluate the product of the Jacobian with a 
// vector v. If NULL, a difference quotient approximation will be used.
void newton_solver_set_jacobian_vector_product(newton_solver_t* solver,
                                               newton_solver_Jv_func Jv_func);

// Call this function to have the Newton solver take a full Newton step on 
// the next solve.
void newton_solver_use_full_step(newton_solver_t* solver);

// Call this function to have the Newton solver use a line search to determine 
// the Newton step size on the next solve.
void newton_solver_use_line_search(newton_solver_t* solver);

// Call this function to have the Newton solver use a fixed point iteration 
// with Anderson acceleration on the next solve. The given (non-negative) 
// number of residuals will be stored and used in the acceleration.
void newton_solver_use_fixed_point(newton_solver_t* solver, 
                                   int num_residuals);

// Call this function to have the Newton solver use a Picard iteration with 
// Anderson acceleration on the next solve. This requires that a Jacobian-vector
// product be given to form a linear matrix L for the iteration. The given 
// (non-negative) number of residuals will be stored and used in the 
// acceleration.
void newton_solver_use_picard(newton_solver_t* solver, 
                              newton_solver_Jv_func Jv_func,
                              int num_residuals);

// Sets the tolerances for the function norm (norm_tolerance) and the Newton
// step (step_tolerance) for the nonlinear solver. The default values 
// for these tolerances are the machine roundoff threshold raised to the 
// 1/3 and 1/2 powers, respectively.
void newton_solver_set_tolerances(newton_solver_t* solver, 
                                  real_t norm_tolerance, 
                                  real_t step_tolerance);

// Sets the maximum number of Newton iterations for the solver.
void newton_solver_set_max_iterations(newton_solver_t* solver, int max_iterations);

// Sets the stopping critieria for the underlying linear solver.
// The options are NEWTON_EISENSTAT_WALKER1, NEWTON_EISENSTAT_WALKER2, and 
// NEWTON_CONSTANT_ETA. See the Kinsol documentation for details on these criteria.
// The default setting is NEWTON_EISENSTAT_WALKER1.
typedef enum
{
  NEWTON_EISENSTAT_WALKER1,
  NEWTON_EISENSTAT_WALKER2,
  NEWTON_CONSTANT_ETA,
} newton_solver_stopping_criteria_t;
void newton_solver_set_linear_solver_stopping_criteria(newton_solver_t* solver,
                                                       newton_solver_stopping_criteria_t criteria);

// In the case that the linear solver is set to NEWTON_CONSTANT_ETA, this function 
// sets the value of the constant coefficient used to determine convergence of the 
// linear solver. See Kinsol documentation for details.
void newton_solver_set_constant_eta(newton_solver_t* solver, real_t eta);

// Gets an internal pointer to the preconditioner.
newton_pc_t* newton_solver_preconditioner(newton_solver_t* solver);

// Evaluates the residual vector, storing it in F.
void newton_solver_eval_residual(newton_solver_t* solver, real_t t, real_t* X, real_t* F);

// Solves the nonlinear system of equations F(X, t) = 0 in place, 
// using X as the initial guess. Returns true if the solution was obtained, 
// false if not. The number of nonlinear iterations will be stored in 
// num_iterations upon success.
bool newton_solver_solve(newton_solver_t* solver, real_t t, real_t* X, int* num_iterations);

// Resets the internal state of the solver to its original state at time t.
void newton_solver_reset(newton_solver_t* solver, real_t t);

// Diagnostics for the nonlinear solver.
typedef struct
{
  char* status_message; // borrowed pointer from solver: do not free.
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
} newton_solver_diagnostics_t;

// Retrieve diagnostics for the nonlinear solver.
void newton_solver_get_diagnostics(newton_solver_t* solver, 
                                   newton_solver_diagnostics_t* diagnostics);

// Writes nonlinear solver diagnostics to the given file.
void newton_solver_diagnostics_fprintf(newton_solver_diagnostics_t* diagnostics, 
                                       FILE* stream);

#endif

