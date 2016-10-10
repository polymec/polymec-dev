// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ARK_ODE_INTEGRATOR_H
#define POLYMEC_ARK_ODE_INTEGRATOR_H

#include "core/krylov_solver.h"
#include "integrators/ode_integrator.h"
#include "integrators/newton_pc.h"

// This type of ODE integrator integrates a multi-scale set of ordinary 
// differential equations using semi-implicit Additive Runge Kutte (ARK) 
// methods. The equations integrated take the form
//
// dU/dt = fe(t, U) + fi(t, U)
//
// where fe is a slowly-varying function that can be integrated explicitly, 
// and fi is a quickly-varying function that should be integrated implicitly.

// This type identifies a predictor method to use for implicit integration.
typedef enum
{
  ARK_TRIVIAL_PREDICTOR, // "trivial predictor" -- uses previous time solution
  ARK_MAXORDER_PREDICTOR, // predictor with maximum order polynomial interpolant
  ARK_VARORDER_PREDICTOR, // predictor with variable order polynomial interpolant
  ARK_CUTOFF_PREDICTOR, // predictor with "cutoff" order -- max order in first half of dt,
                        //                                  first order in last half
  ARK_BOOTSTRAP_PREDICTOR // predictor that uses only current step information
} ark_predictor_t;

// This type identifies a Krylov method to use for the Jacobian-Free
// Newton-Krylov version of this integrator.
typedef enum 
{
  JFNK_ARK_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_ARK_FGMRES,   // Flexible GMRES Krylov solver
  JFNK_ARK_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_ARK_TFQMR,    // Transpose-Free QMR Krylov solver
  JFNK_ARK_PCG       // Preconditioned Conjugate Gradient Krylov solver
} jfnk_ark_krylov_t;

// This constructs an integrator for the slowly-varying system
// dx/dt = fe(t, U) with the desired order of time accuracy.
// The optional stable_dt_func argument supplies a function that computes the 
// next timestep subject to stability constraints.
ode_integrator_t* explicit_ark_ode_integrator_new(int order, 
                                                  MPI_Comm comm,
                                                  int num_local_values, 
                                                  int num_remote_values, 
                                                  void* context, 
                                                  int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                                  real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                                  void (*dtor)(void* context));

// This constructs an integrator for the system with given fe, fi,
// with the desired order of time accuracy, using a fixed-point (functional) 
// iteration method that does not require a Newton-Krylov solver. At least 
// one of fe and fi must be non-NULL. If fi is NULL, the system is assumed to 
// vary "slowly" and will be explicitly integrated; if fe is NULL, the system 
// is assumed to be stiff, and will be implicitly integrated; if both are
// non-NULL, they will be integrated using an adaptive IMEX method.
// max_anderson_accel_dim is the maximum dimension of the underlying Anderson 
// acceleration subspace.
ode_integrator_t* functional_ark_ode_integrator_new(int order, 
                                                    MPI_Comm comm,
                                                    int num_local_values, 
                                                    int num_remote_values, 
                                                    void* context, 
                                                    int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                                    int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                                    real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                                    void (*dtor)(void* context),
                                                    int max_anderson_accel_dim);

// This constructs an integrator for the system with given fe, fi,
// with the desired order of time accuracy, using the Jacobian-Free 
// Newton-Krylov solver of the given type. If the Jacobian-vector 
// product Jy is not given for this integrator, this product will be 
// approximated using difference quotients. fi must be non-NULL, and fe can 
// be either non-NULL or NULL. If fe is NULL, the system is assumed to be 
// stiff, and will be implicitly integrated. The fi_is_linear and 
// fi_is_time_dependent flags provide information about fi in order to help 
// the integrator optimize its performance (and fi_is_time_dependent is 
// ignored if fi is not linear in U).
ode_integrator_t* jfnk_ark_ode_integrator_new(int order, 
                                              MPI_Comm comm,
                                              int num_local_values, 
                                              int num_remote_values, 
                                              void* context, 
                                              int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                              int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                              bool fi_is_linear,
                                              bool fi_is_time_dependent,
                                              real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                              int (*Jy_func)(void* context, real_t t, real_t* U, real_t* fi, real_t* y, real_t* temp, real_t* Jy),
                                              void (*dtor)(void* context),
                                              newton_pc_t* precond,
                                              jfnk_ark_krylov_t solver_type,
                                              int max_krylov_dim);

// This function creates an ARK integrator that uses an inexact Newton-Krylov 
// method to solve the underlying linearized equations. This method constructs 
// a full Jacobian matrix and updates it only when needed, and requires a 
// function to be specified for the update. The given Krylov factory is used 
// to create solver, preconditioner, matrix and vector objects for solving 
// the underlying linear system, and the sparsity pattern is used for the 
// creation of the Jacobian matrix. The integrator assumes control of all 
// these objects.
// The function rhs_func computes the time derivative of the solution, much as 
// it does for jfnk_bdf_ode_integrator_new above. The required function J_func 
// computes the elements of the Jacobian matrix using t, U, and U_dot. The 
// matrix must also be assembled within this function.
ode_integrator_t* ink_ark_ode_integrator_new(int order, 
                                             MPI_Comm comm,
                                             krylov_factory_t* factory,
                                             matrix_sparsity_t* J_sparsity,
                                             void* context, 
                                             int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                             int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                             real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                             bool fi_is_linear,
                                             bool fi_is_time_dependent,
                                             int (*J_func)(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J),
                                             void (*dtor)(void* context));

// Specifies that the INK ARK integrator should use the Preconditioned 
// Conjugate Gradient (PCG) method.
void ink_ark_ode_integrator_use_pcg(ode_integrator_t* ink_ark_ode_integ);

// Specifies that the INK ARK integrator should use the Generalized 
// Minimum Residual (GMRES) method with the specified maximum Krylov subspace
// dimension.
void ink_ark_ode_integrator_use_gmres(ode_integrator_t* ink_ark_ode_integ,
                                      int max_krylov_dim);

// Specifies that the INK ARK integrator should use the Stabilized 
// Bi-Conjugate Gradient (BiCGSTAB) method.
void ink_ark_ode_integrator_use_bicgstab(ode_integrator_t* ink_ark_ode_integ);

// Specifies that the INK ARK Integrator should use the given "special" 
// Krylov solver with the given options.
void ink_ark_ode_integrator_use_special(ode_integrator_t* ink_ark_ode_integ,
                                        const char* solver_name,
                                        string_string_unordered_map_t* options);

// Specifies that the INK ARK integrator should use the preconditioner with 
// the given name, set with the given options.
void ink_ark_ode_integrator_set_pc(ode_integrator_t* ink_ark_ode_integ,
                                   const char* pc_name, 
                                   string_string_unordered_map_t* options);

// Sets the block size for the INK ARK integrator.
void ink_ark_ode_integrator_set_block_size(ode_integrator_t* ink_ark_ode_integ,
                                           int block_size);

// Returns the context pointer for the given INK ARK integrator.
void* ink_ark_ode_integrator_context(ode_integrator_t* ink_ark_ode_integ);

// Convergence failure status codes, used by the ARK integrator machinery below to determine whether 
// the Newton operator needs to be recomputed.
typedef enum
{
  ARK_CONV_NO_FAILURES,   // Local error test failed at previous Newton step, but iteration converged
  ARK_CONV_BAD_J_FAILURE, // Previous Newton corrector iteration did not converge OR linear solve failed in 
                          // A recoverable manner (Jacobian needs updating either way)
  ARK_CONV_OTHER_FAILURE, // Newton iteration failed to converge with current Jacobian data
} ark_conv_status_t;

// This function creates an ARK integrator that uses the given methods to define how
// it solves the linear systems that underlie the ARK method. These methods are:
// * rhs_func -- Used to compute the right-hand side of the ODE.
// * reset_func -- Used to reset the state of the integrator to integrate U at time t, as 
//                 invoked by ode_integrator_reset(integ, t, U).
// * setup_func -- Used to calculate and store a representation of the linear operator 
//                 I - gamma * J, where I is the identity operator, J is the Jacobian, and 
//                 gamma is a positive scale factor related to the time step. This operator
//                 is computed whenever the integrator deems it necessary to reflect the 
//                 the state of the system (according to inexact-Newton-like methods).
//                 Arguments are:
//                 - context: A pointer to the context object that stores the state of the integrator.
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
//                 - context: A pointer to the context object that stores the state of the integrator.
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
//                 - num_iters: A pointer that will store the number of iterations needed to solve 
//                              the linear system.
//                 Should return 0 on success, a positive value for a recoverable error, and a negative value 
//                 for an unrecoverable error.
// * dtor -- Used to destroy the context pointer.
ode_integrator_t* ark_ode_integrator_new(const char* name,
                                         int order, 
                                         MPI_Comm comm,
                                         int num_local_values, 
                                         int num_remote_values, 
                                         void* context, 
                                         int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                         int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                         real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                         bool fi_is_linear,
                                         bool fi_is_time_dependent,
                                         int (*reset_func)(void* context, real_t t, real_t* U),
                                         int (*setup_func)(void* context, 
                                                           ark_conv_status_t conv_status, 
                                                           real_t gamma, 
                                                           int step,
                                                           real_t t, 
                                                           real_t* U_pred, 
                                                           real_t* fi_pred, 
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

// This returns the context pointer passed to the ark_ode_integrator 
// constructor. In general, this will NOT return the same pointer as 
// ode_integrator_context!
void* ark_ode_integrator_context(ode_integrator_t* integrator);

// Sets control parameters for adaptive time stepping. Specifically:
// max_growth (> 1) - The maximum factor by which the step is allowed to grow 
//                    between consecutive steps.
// max_initial_growth (> 1) - The maximum factor by which the step is allowed to grow 
//                            after the first step.
// max_convergence_cut_factor (< 1) - The maximum factor by which the time step is 
//                                    modified after successive convergence failures.
// max_accuracy_cut_factor (< 1) - The maximum factor by which the time step is 
//                                 modified after successive accuracy failures.
// safety_factor - The factor by which the accuracy-based time step is multiplied
//                 by the integrator in taking a step.
// cfl_fraction (<= 1) - The factor by which the stability-based time step is multiplied
//                       by the integrator in taking a step (ignored if there's no fe).
void ark_ode_integrator_set_step_controls(ode_integrator_t* integrator,
                                          real_t max_growth,
                                          real_t max_initial_growth,
                                          real_t max_convergence_cut_factor,
                                          real_t max_accuracy_cut_factor,
                                          real_t safety_factor,
                                          real_t cfl_fraction);

// Sets the method to use for predicting implicit terms in the RK integration.
void ark_ode_integrator_set_predictor(ode_integrator_t* integrator, 
                                      ark_predictor_t predictor);

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

// Evaluates the explicit part of the system at the given time and with the 
// given solution X, placing the results in fe.
void ark_ode_integrator_eval_fe(ode_integrator_t* integrator, real_t t, real_t* X, real_t* fe);

// Evaluates the implicit part of the system at the given time and with the 
// given solution X, placing the results in fi.
void ark_ode_integrator_eval_fi(ode_integrator_t* integrator, real_t t, real_t* X, real_t* fi);

// Returns a stable timestep for the solution X at time t.
real_t ark_ode_integrator_stable_dt(ode_integrator_t* integrator, real_t t, real_t* X);

// Returns an internal pointer to the preconditioner passed to this 
// integrator during construction time.
newton_pc_t* ark_ode_integrator_preconditioner(ode_integrator_t* integrator);

// Diagnostics for the time integrator.
typedef struct
{
  char* status_message; // borrowed pointer from integrator: do not free.
  long int num_steps, num_explicit_steps, num_accuracy_steps,
           num_step_attempts;
  real_t initial_step_size, last_step_size, next_step_size;
  real_t t;
  long int num_fe_evaluations, num_fi_evaluations;
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
                                         void (*fe_computed)(void* context, real_t t, real_t* U, real_t* U_dot),
                                         void (*fi_computed)(void* context, real_t t, real_t* U, real_t* U_dot),
                                         void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* Jy),
                                         void (*dtor)(void* context));

// Adds the given observer to the given am_ode_integrator. The observer 
// is consumed by the integrator, so no destructor is needed.
void ark_ode_integrator_add_observer(ode_integrator_t* integrator,
                                     ark_ode_observer_t* observer);

#endif

