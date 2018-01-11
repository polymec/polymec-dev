// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BDF_ODE_SOLVER_H
#define POLYMEC_BDF_ODE_SOLVER_H

#include "solvers/krylov_solver.h"
#include "solvers/ode_solver.h"
#include "solvers/newton_pc.h"

// This type of ODE solver integrates a stiff set of ordinary 
// differential equations using Backwards Difference Formulae (BDF). These 
// formulae can provide accurate time integration at an order from 1 to 5, and 
// L stability at orders 1 and 2. The BDF solvers all use Newton iteration.
// These solvers are implemented using CVODE from the Sundials suite 
// of nonlinear solvers.

// Variants of the Jacobian-Free Newton-Krylov BDF solvers described below.
typedef enum 
{
  JFNK_BDF_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_BDF_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_BDF_TFQMR     // Transpose-Free QMR Krylov solver
} jfnk_bdf_krylov_t;

// This function creates a BDF solver that uses a Jacobian-Free
// Newton-Krylov method to solve the underlying linearized equations. This 
// method requires a preconditioner that is a coarse approximation of the 
// Jacobian matrix, but captures its essential behavior. The right-hand side 
// function and destructor are required, and a function to compute Jy, the 
// product of the Jacobian with the vector y, can be optionally provided. Its 
// signature is: 
// Jy_func(context, t, U, U_dot, y, Jy, temp), where context is the context pointer, 
// t is the time at which Jy is evaluated, y is the vector the Jacobian operator
// J is applied to, U_dot is the right hand side function evaluated at t and U,
// temp is a work vector the same size as y, and Jy is a vector that stores the 
// product Jy.
// If Jy_func is not given, a finite difference approximation of Jy will be used.
// Additionally, the type of Krylov solver (JFNK_BDF_GMRES, JFNK_BDF_BICGSTAB, 
// or JFNK_BDF_TFQMR) must be given, along with the maximum dimension of the 
// Krylov subspace. 
ode_solver_t* jfnk_bdf_ode_solver_new(int order, 
                                      MPI_Comm comm,
                                      int num_local_values, 
                                      int num_remote_values, 
                                      void* context, 
                                      int (*rhs_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                                      int (*Jy_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* Jy, real_t* temp),
                                      void (*dtor)(void* context),
                                      newton_pc_t* precond,
                                      jfnk_bdf_krylov_t solver_type,
                                      int max_krylov_dim);

// This function creates a BDF solver that uses an inexact Newton-Krylov 
// method to solve the underlying linearized equations. This method constructs 
// a full Jacobian matrix and updates it only when needed, and requires a 
// function to be specified for the update. The given Krylov factory is used 
// to create solver, preconditioner, matrix and vector objects for solving 
// the underlying linear system, and the sparsity pattern is used for the 
// creation of the Jacobian matrix. The solver assumes control of all 
// these objects.
// The function rhs_func computes the time derivative of the solution, much as 
// it does for jfnk_bdf_ode_solver_new above. The required function J_func 
// computes the elements of the Jacobian matrix using t, U, and U_dot. The 
// matrix must also be assembled within this function.
ode_solver_t* ink_bdf_ode_solver_new(int order, 
                                     MPI_Comm comm,
                                     krylov_factory_t* factory,
                                     matrix_sparsity_t* J_sparsity,
                                     void* context, 
                                     int (*rhs_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                                     int (*J_func)(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J),
                                     void (*dtor)(void* context));

// Specifies that the INK BDF solver should use the Preconditioned 
// Conjugate Gradient (PCG) method.
void ink_bdf_ode_solver_use_pcg(ode_solver_t* ink_bdf_ode_integ);

// Specifies that the INK BDF solver should use the Generalized 
// Minimum Residual (GMRES) method with the specified maximum Krylov subspace
// dimension.
void ink_bdf_ode_solver_use_gmres(ode_solver_t* ink_bdf_ode_integ,
                                  int max_krylov_dim);

// Specifies that the INK BDF solver should use the Stabilized 
// Bi-Conjugate Gradient (BiCGSTAB) method.
void ink_bdf_ode_solver_use_bicgstab(ode_solver_t* ink_bdf_ode_integ);

// Specifies that the INK BDF solver should use the given "special" 
// Krylov solver with the given options.
void ink_bdf_ode_solver_use_special(ode_solver_t* ink_bdf_ode_integ,
                                    const char* solver_name,
                                    string_string_unordered_map_t* options);

// Specifies that the INK BDF solver should use the preconditioner with 
// the given name, set with the given options.
void ink_bdf_ode_solver_set_pc(ode_solver_t* ink_bdf_ode_integ,
                               const char* pc_name, 
                               string_string_unordered_map_t* options);

// Sets the block size for the INK BDF solver.
void ink_bdf_ode_solver_set_block_size(ode_solver_t* ink_bdf_ode_integ,
                                       int block_size);

// Returns the context pointer for the given INK BDF solver.
void* ink_bdf_ode_solver_context(ode_solver_t* ink_bdf_ode_integ);

// Convergence failure status codes, used by the BDF solver machinery below to determine whether 
// the Newton operator needs to be recomputed.
typedef enum
{
  BDF_CONV_NO_FAILURES,   // Local error test failed at previous Newton step, but iteration converged
  BDF_CONV_BAD_J_FAILURE, // Previous Newton corrector iteration did not converge OR linear solve failed in 
                          // A recoverable manner (Jacobian needs updating either way)
  BDF_CONV_OTHER_FAILURE, // Newton iteration failed to converge with current Jacobian data
} bdf_conv_status_t;

// This function creates a BDF solver that uses the given methods to define how
// it solves the linear systems that underlie the BDF method. These methods are:
// * rhs_func -- Used to compute the right-hand side of the ODE.
// * reset_func -- Used to reset the state of the solver to integrate U at time t, as 
//                 invoked by ode_solver_reset(integ, t, U).
// * setup_func -- Used to calculate and store a representation of the linear operator 
//                 I - gamma * J, where I is the identity operator, J is the Jacobian, and 
//                 gamma is a positive scale factor related to the time step. This operator
//                 is computed whenever the solver deems it necessary to reflect the 
//                 the state of the system (according to inexact-Newton-like methods).
//                 Arguments are:
//                 - context: A pointer to the context object that stores the state of the solver.
//                 - conv_status: A status code produced by the Newton solver for its solution
//                                at the given time step. See above.
//                 - gamma: the scaling factor in I - gamma * J.
//                 - step: The current integration step number since the initial time.
//                 - t: The current time.
//                 - U_pred: The predicted solution vector for the current step.
//                 - U_dot_pred: The value of the right hand side at time t and U = U_pred.
//                 - J_updated: A pointer to a boolean variable that conveys whether J has been updated by 
//                              the function.
//                 - work1, work2, work3: work vectors of the same size as the solution vector, provided 
//                                        for use by this method.
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * solve_func -- Used to solve the linear system (I - gamma * J) * X = B. Arguments are:
//                 - context: A pointer to the context object that stores the state of the solver.
//                 - t: The current time.
//                 - U: The current solution at time t.
//                 - U_dot: The current right hand side vector for the ODE at time t.
//                 - W: A vector containing error weights, which can be used to enable to computation of 
//                      weighted norms used to test for convergence of any iterative methods within 
//                      the solver.
//                 - res_norm_tol: The tolerance to set on the residual norm to consider the system 
//                                 successfully solved.
//                 - B: The right hand side vector for the linear system, which will be replaced by 
//                      X, the solution to the linear system.
//                 - num_iters: A pointer that will store the number of iterations taken to solve 
//                              the linear system.
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * dtor -- Used to destroy the context pointer.
ode_solver_t* bdf_ode_solver_new(const char* name,
                                 int order, 
                                 MPI_Comm comm,
                                 int num_local_values, 
                                 int num_remote_values, 
                                 void* context, 
                                 int (*rhs_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                                 int (*reset_func)(void* context, real_t t, real_t* U),
                                 int (*setup_func)(void* context, 
                                                   bdf_conv_status_t conv_status, 
                                                   real_t gamma, 
                                                   int step,
                                                   real_t t, 
                                                   real_t* U_pred, 
                                                   real_t* U_dot_pred, 
                                                   bool* J_updated, 
                                                   real_t* work1, real_t* work2, real_t* work3),
                                 int (*solve_func)(void* context, 
                                                   real_t t, 
                                                   real_t* U,
                                                   real_t* U_dot,
                                                   real_t* W, 
                                                   real_t res_norm_tol,
                                                   real_t* B, 
                                                   int* num_iters), 
                                 void (*dtor)(void* context));

// This returns the context pointer passed to the bdf_ode_solver 
// constructor that created this solver. In general, this will NOT 
// return the same pointer as ode_solver_context!
void* bdf_ode_solver_context(ode_solver_t* solver);

// Sets the relative and absolute tolerances for integrated quantities.
void bdf_ode_solver_set_tolerances(ode_solver_t* solver,
                                   real_t relative_tol, real_t absolute_tol);

// Sets the error weights for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances. Weights are copied into the solver.
void bdf_ode_solver_set_error_weights(ode_solver_t* solver, real_t* weights);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void bdf_ode_solver_set_error_weight_function(ode_solver_t* solver,
                                              void (*compute_weights)(void* context, real_t* y, real_t* weights));

// Toggles stability limit detection, making the BDF method more robust at 
// a cost of 10-20% overhead.
void bdf_ode_solver_set_stability_limit_detection(ode_solver_t* solver,
                                                  bool use_detection);

// Sets the maximum number of error test failures permitted in attempting 
// a single time step. By default, this value is 7.
void bdf_ode_solver_set_max_err_test_failures(ode_solver_t* solver,
                                              int max_failures);

// Sets the maximum number of nonlinear solver iterations per time step.
// By default, this value is 3.
void bdf_ode_solver_set_max_nonlinear_iterations(ode_solver_t* solver,
                                                 int max_iterations);

// Sets the safety factor (coefficient) used in the nonlinear convergence test.
// By default, this value is 0.1.
void bdf_ode_solver_set_nonlinear_convergence_coeff(ode_solver_t* solver,
                                                    real_t coefficient);

// Evaluates the right-hand side of the system at the given time and with the 
// given solution U, placing the results in U_dot.
void bdf_ode_solver_eval_rhs(ode_solver_t* integ, real_t t, real_t* U, real_t* U_dot);

// Returns an internal pointer to the preconditioner passed to this 
// solver during construction time.
newton_pc_t* bdf_ode_solver_preconditioner(ode_solver_t* solver);

// Diagnostics for the time solver.
typedef struct
{
  char* status_message; // borrowed pointer from solver: do not free.
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
} bdf_ode_solver_diagnostics_t;

// Retrieve diagnostics for the time solver.
void bdf_ode_solver_get_diagnostics(ode_solver_t* solver, 
                                    bdf_ode_solver_diagnostics_t* diagnostics);

// Writes time solver diagnostics to the given file.
void bdf_ode_solver_diagnostics_fprintf(bdf_ode_solver_diagnostics_t* diagnostics, 
                                        FILE* stream);

// This observer type can be used to define objects that respond to actions
// taken by the bdf_ode_solver.
typedef struct bdf_ode_observer_t bdf_ode_observer_t;

// Creates and returns a newly-allocated observer that observes an 
// am_ode_solver. Responses are given by the following arguments:
// rhs_computed - This function is called when the right hand side of the ODE
//                system is computed by the solver.
// Jy_computed - This function is called when the Jacobian-vector product J*y 
//               is computed by the solver.
// Both of these functions are fed the given context object.
bdf_ode_observer_t* bdf_ode_observer_new(void* context,
                                         void (*rhs_computed)(void* context, real_t t, real_t* U, real_t* U_dot),
                                         void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* Jy),
                                         void (*dtor)(void* context));

// Adds the given observer to the given am_ode_solver. The observer 
// is consumed by the solver, so no destructor is needed.
void bdf_ode_solver_add_observer(ode_solver_t* solver,
                                 bdf_ode_observer_t* observer);

#endif

