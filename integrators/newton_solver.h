// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NEWTON_SOLVER_H
#define POLYMEC_NEWTON_SOLVER_H

#include "core/polymec.h"
#include "core/krylov_solver.h"
#include "integrators/newton_pc.h"

// This class represents a collection of algorithms for integrating partial 
// differential equations that are discretized into a (sparse) system of 
// nonlinear equations. The integration is performed using Matrix-free 
// Newton-Krylov methods provided by KINSol with right preconditioning.
typedef struct newton_solver_t newton_solver_t;

// This type distinguishes between the different strategies used to perform
// the Newton iteration.
typedef enum
{
  NEWTON_FULL_STEP,   // Full Newton step
  NEWTON_LINE_SEARCH, // Line search
  NEWTON_PICARD,      // Picard iteration
  NEWTON_FP           // Fixed-point iteration
} newton_solver_strategy_t; 

// Creates a generic Newton solver that uses the given methods to define how it 
// solves the linear systems that underlie the Newton iteration. These methods 
// are:
// * F_func -- Used to compute the system function F(t, U).
// * reset_func -- Used to reset the state of the solver to a different state.
// * setup_func -- Used to calculate and store a representation of the Jacobian matrix
//                 J. This operator is computed whenever the integrator deems it necessary to 
//                 reflect the the state of the system (according to inexact-Newton-like methods).
//                 Arguments are:
//                 - context: A pointer to the context object that stores the state of the solver.
//                 - strategy: The type of iteration being used to solve the system. This affects
//                             what is meant by "J": In the case of full Newton step and line 
//                             search, J = dF/dU. In the case of Picard iteration, J means the 
//                             linearization L of J. This function is never called in the case
//                             of fixed point iteration.
//                 - new_U: A flag that indicates whether U has been updated 
//                          since the last call to the set-up function.
//                 - t: The current time.
//                 - U: The current solution vector U.
//                 - F: The current system function (residual) vector F.
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * solve_func -- Used to solve the linear system J * p = -F. Arguments are:
//                 - context: A pointer to the context object that stores the state of the solver.
//                 - DF: A vector containing the scaling vector DF.
//                 - t: The current time.
//                 - B: The right hand side vector B = -F for the system at time t.
//                 - res_norm_tol: The tolerance against which the residual 2-norm ||-F-J*p||_2 
//                                 will be measured.
//                 - p: An array containing the solution to the linear system, which is the 
//                      increment p = (U(n) - U(n-1)). On input it contains an initial guess 
//                      for the increment; on output it contains the actual increment.
//                 - Jp_norm: A pointer that will store the scaled L2 norm of the product Jp:
//                            ||DF*J*p||_2.
//                 - F_o_Jp: A pointer that will store the scaled dot product of F with Jp:
//                           (DF*F) o (DF*J*p).
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * dtor -- Used to destroy the context pointer.
newton_solver_t* newton_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                   int (*reset_func)(void* context),
                                   int (*setup_func)(void* context, 
                                                     newton_solver_strategy_t strategy,
                                                     bool new_U,
                                                     real_t t,
                                                     real_t* U,
                                                     real_t* F),
                                   int (*solve_func)(void* context, 
                                                     real_t* DF, 
                                                     real_t t, 
                                                     real_t* B,
                                                     real_t res_norm_tol,
                                                     real_t* p,
                                                     real_t* Jp_norm, 
                                                     real_t* F_o_Jp), 
                                   void (*dtor)(void* context));

// This type represents the different Krylov methods that can be used in 
// Jacobian-Free Newton-Krylov (JFNK) implementations of Newton solvers.
typedef enum 
{
  NEWTON_GMRES,    // Generalized minimum residual Krylov solver
  NEWTON_FGMRES,   // Flexible GMRES Krylov solver
  NEWTON_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  NEWTON_TFQMR     // Transpose-Free QMR Krylov solver
} jfnk_newton_t;

// Creates a Newton solver that uses a Jacobian-Free Newton-Krylov (JFNK) 
// method with a given maximum subspace dimension of max_krylov_dim. If the 
// solver_type is NEWTON_GMRES or NEWTON_FGMRES, the maximum number of restarts 
// is given by max_restarts--otherwise that parameter is ignored. 
// The function F_func computes the system function F(U) = R. The function Jv_func, 
// if provided, computes the product of the Jacobian matrix J with a vector v at time t, 
// given a solution U, and the flag new_U set to true if the solution has been updated 
// since the last call (and false if not).
newton_solver_t* jfnk_newton_solver_new(MPI_Comm comm,
                                        int num_local_values,
                                        int num_remote_values,
                                        void* context,
                                        int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                        int (*Jv_func)(void* context, bool new_U, real_t t, real_t* U, real_t* v, real_t* Jv),
                                        void (*dtor)(void* context),
                                        newton_pc_t* precond,
                                        jfnk_newton_t solver_type,
                                        int max_krylov_dim,
                                        int max_restarts);

// This function creates a Newton solver that uses an inexact Newton-Krylov 
// method to solve the underlying linearized equations. This method constructs 
// a full matrix and updates it only when needed, and requires a function to be 
// specified for the update. The given Krylov factory is used to create solver, 
// preconditioner, matrix and vector objects for solving the underlying linear 
// system, and the sparsity pattern is used for the creation of the matrix. The 
// integrator assumes control of all these objects.
// The function F_func computes the system function, much as it does for 
// jfnk_newton_solver_new above. The required function J_func computes the 
// elements of the Jacobian/Picard matrix using t and U. The matrix must also 
// be assembled within this function.
newton_solver_t* ink_newton_solver_new(MPI_Comm comm,
                                       krylov_factory_t* factory,
                                       matrix_sparsity_t* J_sparsity,
                                       void* context, 
                                       int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                       int (*J_func)(void* context, real_t t, real_t* U, real_t* F, krylov_matrix_t* J),
                                       void (*dtor)(void* context));

// Frees a solver.
void newton_solver_free(newton_solver_t* solver);

// Returns the context pointer for the solver.
void* newton_solver_context(newton_solver_t* solver);

// Returns the number of (local) equations in the nonlinear system.
int newton_solver_num_equations(newton_solver_t* solver);

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
// Anderson acceleration on the next solve. The given (non-negative) number 
// of residuals will be stored and used in the acceleration.
void newton_solver_use_picard(newton_solver_t* solver, 
                              int num_residuals);

// Sets the scaling coefficients for the solution vector U such that DU * U 
// all have roughly the same magnitude as F(x) approaches 0, where DU is 
// a diagonal matrix whose elements are the scaling coefficients. If DU is 
// NULL, no scaling is performed.
void newton_solver_set_U_scale(newton_solver_t* solver, 
                               real_t* DU);

// Sets the scaling coefficients for the system function vector F such that 
// DF * F(t, U) all have roughly the same magnitude as F(x) approaches 0, 
// where DF is a diagonal matrix whose elements are the scaling coefficients.
// If DF is NULL, no scaling is performed.
void newton_solver_set_F_scale(newton_solver_t* solver, 
                               real_t* DF);

// Sets algebraic constraints on the components of the solution vector U. 
// If constraints[i] is:
//   0.0  - no constraint is placed on U[i].
//   1.0  - U[i] must be non-negative.
//  -1.0  - U[i] must be non-positive.
//   2.0  - U[i] must be positive.
void newton_solver_set_constraints(newton_solver_t* solver,
                                   real_t* constraints);

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

// Evaluates the residual vector, storing it in R.
void newton_solver_eval_residual(newton_solver_t* solver, 
                                 real_t t, 
                                 real_t* U, 
                                 real_t* R);

// Solves the nonlinear system of equations F(U, t) = 0 in place, 
// using U as the initial guess. Returns true if the solution was obtained, 
// false if not. The number of nonlinear iterations will be stored in 
// num_iterations upon success.
bool newton_solver_solve(newton_solver_t* solver, 
                         real_t t, 
                         real_t* U, 
                         int* num_iterations);

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

