// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DAE_SOLVER_H
#define POLYMEC_DAE_SOLVER_H

#include "solvers/krylov_solver.h"
#include "solvers/newton_pc.h"

// This type indicates whether a given equation in a DAE system is 
// algebraic or differential.
typedef enum
{
  DAE_ALGEBRAIC,
  DAE_DIFFERENTIAL
} dae_equation_t;

// This can be passed as an argument to dae_solver_new to indicate that 
// all equations are algebraic.
extern dae_equation_t* DAE_ALL_ALGEBRAIC;

// This can be passed as an argument to dae_solver_new to indicate that 
// all equations are differential.
extern dae_equation_t* DAE_ALL_DIFFERENTIAL;

// This type indicates constraints (if any) on equations in a DAE system.
typedef enum
{
  DAE_UNCONSTRAINED,
  DAE_NEGATIVE,
  DAE_NONPOSITIVE,
  DAE_NONNEGATIVE,
  DAE_POSITIVE
} dae_constraint_t;

// This can be passed as an argument to dae_solver_new to indicate that 
// no equations have constraints.
extern dae_constraint_t* DAE_ALL_UNCONSTRAINED;

// This can be passed as an argument to dae_solver_new to indicate that 
// all equations have solution values that are negative.
extern dae_constraint_t* DAE_ALL_NEGATIVE;

// This can be passed as an argument to dae_solver_new to indicate that 
// all equations have solution values that are non-positive.
extern dae_constraint_t* DAE_ALL_NONPOSITIVE;

// This can be passed as an argument to dae_solver_new to indicate that 
// all equations have solution values that are non-negative.
extern dae_constraint_t* DAE_ALL_NONNEGATIVE;

// This can be passed as an argument to dae_solver_new to indicate that 
// all equations have solution values that are positive.
extern dae_constraint_t* DAE_ALL_POSITIVE;

// This type describes methods of correcting initial conditions (U, U_dot) for 
// differential-algebraic systems.
typedef enum
{
  DAE_IC_CORRECT_DERIVATIVES, // Compute corrections to U_dot given U.
  DAE_IC_ASSUME_QUASISTATIC,  // Compute corrections to U assuming U_dot = 0.
  DAE_IC_ASSUME_CONSISTENT    // Assume that U and U_dot are consistent.
} dae_ic_correction_t;

// Types of Krylov solver to use for the Jacobian-Free Newton-Krylov DAE method.
typedef enum 
{
  JFNK_DAE_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_DAE_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_DAE_TFQMR     // Transpose-Free QMR Krylov solver
} jfnk_dae_krylov_t;

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. 
typedef struct dae_solver_t dae_solver_t;

// Creates an solver that uses a Jacobian-Free Newton-Krylov method to 
// solve a system of Differential Algebraic Equations (DAE) with a given maximum 
// subspace dimension of max_krylov_dim. The equation_types array (of length 
// num_local_values) indicates whether each equation in the system is algebraic 
// (DAE_ALGEBRAIC) or differential (DAE_DIFFERENTIAL). The constraints array 
// indicates any constraints on the solution for a given component 
// (DAE_UNCONSTRAINED, DAE_NEGATIVE, DAE_NONPOSITIVE, DAE_NONNEGATIVE, 
// DAE_POSITIVE). See the DAE_ALL* symbols above for shorthand ways of 
// expressing systems with simple equation types and constraints. This method 
// requires a preconditioner that is a coarse approximation of the 
// Newton matrix, but captures its essential behavior. The right-hand side 
// function and destructor are required, and a function to compute Jy, the 
// product of the Jacobian J = dF/dU + alpha * dF/d(U_dot) with the vector y, 
// can be optionally provided. Its signature is: 
// Jy_func(context, t, U, alpha, U_dot, F, y, Jy, tmp1, tmp2), where context is the 
// context pointer, t is the time at which Jy is evaluated, y is the vector the Jacobian 
// operator J is applied to, F is the residual function evaluated at (t, U, U_dot),
// Jy is a vector that stores the product Jy, and tmp1, tmp2 are work vectors the same size 
// as y. If Jy_func is not given, a finite difference approximation of it will be used.
// Additionally, the type of Krylov solver (JFNK_DAE_GMRES, JFNK_DAE_BICGSTAB, 
// or JFNK_DAE_TFQMR) must be given, along with the maximum dimension of the 
// Krylov subspace. 
dae_solver_t* jfnk_dae_solver_new(int order,
                                  MPI_Comm comm,
                                  dae_equation_t* equation_types,
                                  dae_constraint_t* constraints,
                                  int num_local_values,
                                  int num_remote_values,
                                  void* context,
                                  int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F),
                                  int (*Jy_func)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F,
                                                 real_t* y, real_t* Jy, real_t* tmp1, real_t* tmp2),
                                  void (*dtor)(void* context),
                                  newton_pc_t* precond,
                                  jfnk_dae_krylov_t solver_type,
                                  int max_krylov_dim);

// This function creates a DAE solver that uses the given methods to define how
// it solves the linear systems that underlie the BDF method. These methods are:
// * F_func -- Used to compute the system/residual function F(t, U, U_dot).
// * reset_func -- Used to reset the state of the solver to integrate U at time t, as 
//                 invoked by dae_solver_reset(integ, t, U, U_dot).
// * setup_func -- Used to calculate and store a representation of the linear operator 
//                 dF/dU + alpha * dF/d(U_dot), where alpha is a positive scale factor 
//                 related to the time step. This operator is computed whenever the 
//                 solver deems it necessary to reflect the the state of the system 
//                 (according to inexact-Newton-like methods).
//                 Arguments are:
//                 - context: A pointer to the context object that stores the state of the solver.
//                 - alpha: the scaling factor in dF/dU + alpha * dF/d(U_dot).
//                 - step: The current integration step number since the initial time.
//                 - t: The current time.
//                 - U_pred: The predicted solution vector U for the current step.
//                 - U_dot_pred: The predicted time derivative U_dot for the current step.
//                 - F_pred: The residual F(t, U_pred, U_dot_pred).
//                 - work1, work2, work3: work vectors of the same size as the solution vector, provided 
//                                        for use by this method.
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * solve_func -- Used to solve the linear system J * X = B. Arguments are:
//                 - context: A pointer to the context object that stores the state of the solver.
//                 - t: The current time.
//                 - U: The solution at time t.
//                 - U_dot: The solution time derivative the DAE system at time t.
//                 - F: The system function/residual at time t.
//                 - W: A vector containing error weights, which can be used to enable to computation of 
//                      weighted norms used to test for convergence of any iterative methods within 
//                      the solver.
//                 - res_norm_tol: The tolerance to set on the residual norm to consider the system 
//                                 successfully solved.
//                 - B: The right hand side vector for the linear system, which will be replaced by 
//                      X, the solution to the linear system.
//                 - num_iters: A pointer that will store the number of iterations needed to solve
//                              the linear system.
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * dtor -- Used to destroy the context pointer.
dae_solver_t* dae_solver_new(const char* name,
                             int order, 
                             MPI_Comm comm,
                             dae_equation_t* equation_types,
                             dae_constraint_t* constraints,
                             int num_local_values, 
                             int num_remote_values, 
                             void* context, 
                             int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F),
                             int (*reset_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                             int (*setup_func)(void* context, 
                                               real_t alpha, 
                                               int step,
                                               real_t t, 
                                               real_t* U_pred, 
                                               real_t* U_dot_pred, 
                                               real_t* F_pred, 
                                               real_t* work1, real_t* work2, real_t* work3),
                             int (*solve_func)(void* context, 
                                               real_t t, 
                                               real_t* U,
                                               real_t* U_dot,
                                               real_t* F,
                                               real_t* W, 
                                               real_t res_norm_tol,
                                               real_t* B,
                                               int* num_iters), 
                             void (*dtor)(void* context));

// Frees a DAE solver.
void dae_solver_free(dae_solver_t* solver);

// Returns an internal pointer to the string containing the name 
// of this solver.
char* dae_solver_name(dae_solver_t* solver);

// Returns the context object for this solver.
void* dae_solver_context(dae_solver_t* solver);

// Returns the order of the integration method.
int dae_solver_order(dae_solver_t* solver);

// Gets an internal pointer to the preconditioner.
newton_pc_t* dae_solver_preconditioner(dae_solver_t* solver);

// Sets whether to use a stability limit detection algorithm to improve 
// robustness on particularly stiff problems (2-10% overhead, depending 
// on the problem).
void dae_solver_set_stability_limit_detection(dae_solver_t* solver,
                                              bool use_detection);

// Sets the relative and absolute tolerances for integrated quantities.
void dae_solver_set_tolerances(dae_solver_t* solver,
                               real_t relative_tol, real_t absolute_tol);

// Sets the error weights for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances. Weights are copied into the solver.
void dae_solver_set_error_weights(dae_solver_t* solver, real_t* weights);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void dae_solver_set_error_weight_function(dae_solver_t* solver,
                                          void (*compute_weights)(void* context, real_t* y, real_t* weights));

// Evaluates the right-hand side of the system at the given time and with the 
// given solution U, placing the results in rhs.
void dae_solver_eval_rhs(dae_solver_t* integ, real_t t, real_t* U, real_t* rhs);

// Sets the maximum time step size for the next integration step.
void dae_solver_set_max_dt(dae_solver_t* integ, real_t max_dt);

// Sets the time past which the solver will not step.
void dae_solver_set_stop_time(dae_solver_t* integ, real_t stop_time);

// Integrates the given solution U in place (given its time derivative U_dot), 
// taking a single step starting at time *t and storing the new time in *t as 
// well. Also evolves U_dot. Returns true if the step succeeded, false if it 
// failed for some reason. If a step fails, t, U, and U_dot remain unchanged.
bool dae_solver_step(dae_solver_t* solver, real_t max_dt, real_t* t, real_t* U, real_t* U_dot);

// Resets the solver to prepare it to take a step when U, U_dot, and/or t 
// have changed by some process outside of the solver. This resets any 
// history information stored within the solver. Values of U are corrected
// in order to be made consistent with the given values of U_dot.
// If correct_initial_conditions is true, the solver will attempt to 
// correct U given U_dot.
void dae_solver_reset(dae_solver_t* solver, 
                      real_t t, real_t* U, real_t* U_dot,
                      dae_ic_correction_t ic_correction);

// Diagnostics for the time solver.
typedef struct
{
  char* status_message; // borrowed pointer from solver: do not free.
  long int num_steps;
  int order_of_last_step;
  real_t initial_step_size;
  real_t last_step_size;
  long int num_residual_evaluations;
  long int num_linear_solve_setups;
  long int num_linear_solve_iterations;
  long int num_linear_solve_convergence_failures;
  long int num_error_test_failures;
  long int num_nonlinear_solve_iterations;
  long int num_nonlinear_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
} dae_solver_diagnostics_t;

// Retrieve diagnostics for the time solver.
void dae_solver_get_diagnostics(dae_solver_t* solver, 
                                dae_solver_diagnostics_t* diagnostics);

// Writes time solver diagnostics to the given file.
void dae_solver_diagnostics_fprintf(dae_solver_diagnostics_t* diagnostics, 
                                    FILE* stream);

// This function creates a DAE solver that uses an inexact Newton-Krylov 
// method to solve the underlying linearized equations. This method constructs 
// a full Jacobian matrix and updates it only when needed, and requires a 
// function to be specified for the update. The given Krylov factory is used 
// to create solver, preconditioner, matrix and vector objects for solving 
// the underlying linear system, and the sparsity pattern is used for the 
// creation of the Jacobian matrix. The solver assumes control of all 
// these objects.
// The function F_func computes the system function/residual F(t, U, U_dot), much as 
// it does for jfnk_dae_solver_new above. The required function J_func 
// computes the elements of the Jacobian matrix J = dF/dU + alpha * dF/d(U_dot) 
// using t, U, alpha, and U_dot. The matrix must also be assembled within this function.
dae_solver_t* ink_dae_solver_new(int order, 
                                 MPI_Comm comm,
                                 dae_equation_t* equation_types,
                                 dae_constraint_t* constraints,
                                 krylov_factory_t* factory,
                                 matrix_sparsity_t* J_sparsity,
                                 void* context, 
                                 int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F),
                                 int (*J_func)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F, krylov_matrix_t* J),
                                 void (*dtor)(void* context));

// Specifies that the INK DAE solver should use the Preconditioned 
// Conjugate Gradient (PCG) method.
void ink_dae_solver_use_pcg(dae_solver_t* ink_dae_integ);

// Specifies that the INK DAE solver should use the Generalized 
// Minimum Residual (GMRES) method with the specified maximum Krylov subspace
// dimension.
void ink_dae_solver_use_gmres(dae_solver_t* ink_dae_integ,
                              int max_krylov_dim);

// Specifies that the INK DAE solver should use the Stabilized 
// Bi-Conjugate Gradient (BiCGSTAB) method.
void ink_dae_solver_use_bicgstab(dae_solver_t* ink_dae_integ);

// Specifies that the INK DAE solver should use the given "special" 
// Krylov solver with the given options.
void ink_dae_solver_use_special(dae_solver_t* ink_dae_integ,
                                const char* solver_name,
                                string_string_unordered_map_t* options);

// Specifies that the INK DAE solver should use the preconditioner with 
// the given name, set with the given options.
void ink_dae_solver_set_pc(dae_solver_t* ink_dae_integ,
                           const char* pc_name, 
                           string_string_unordered_map_t* options);

// Sets the block size for the INK DAE solver.
void ink_dae_solver_set_block_size(dae_solver_t* ink_dae_integ,
                                   int block_size);

// Returns the context pointer for the given INK DAE solver.
void* ink_dae_solver_context(dae_solver_t* ink_dae_integ);

#endif

