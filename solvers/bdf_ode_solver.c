// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "core/array.h"
#include "core/timer.h"
#include "solvers/bdf_ode_solver.h"
#include "solvers/newton_pc.h"
#include "cvode/cvode.h"

// JFNK stuff.
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunlinsol/sunlinsol_spfgmr.h"
#include "sunlinsol/sunlinsol_spbcgs.h"
#include "sunlinsol/sunlinsol_sptfqmr.h"

#include "sunnonlinsol/sunnonlinsol_newton.h"

// Stuff for generalized BDF solvers.
#include "cvode/cvode_impl.h"

struct bdf_ode_observer_t 
{
  void* context;
  void (*rhs_computed)(void* context, real_t t, real_t* U, real_t* U_dot);
  void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* Jy);
  void (*dtor)(void* context);
};

typedef struct
{
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context; 
  int (*rhs)(void* context, real_t t, real_t* U, real_t* xdot);
  void (*dtor)(void* context);

  // CVODE data structures.
  void* cvode;
  N_Vector U; 
  real_t* U_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.
  SUNLinearSolver ls;
  SUNNonlinearSolver nls;

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  int (*Jy)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* temp, real_t* Jy);

  // Generalized adaptor stuff.
  real_t sqrtN;
  int (*reset_func)(void* context, real_t t, real_t* U);
  int (*setup_func)(void* context, 
                    bdf_conv_status_t conv_status, 
                    real_t gamma, 
                    int step,
                    real_t t, 
                    real_t* U_pred, 
                    real_t* U_dot_pred, 
                    bool* J_current, 
                    real_t* work1, real_t* work2, real_t* work3);
  int (*solve_func)(void* context, 
                    real_t t, 
                    real_t* U,
                    real_t* U_dot,
                    real_t* W, 
                    real_t res_norm_tol,
                    real_t* B,
                    int* num_iters);
  int num_linear_iterations, num_linear_conv_failures;

  // Error weight function.
  void (*compute_weights)(void* context, real_t* y, real_t* weights);
  real_t* error_weights;

  // Observers.
  ptr_array_t* observers;
} bdf_ode_t;

static char* get_status_message(int status, real_t current_time)
{
  char* status_message = NULL;
  if (status == CV_TOO_CLOSE)
    status_message = string_dup("t1 and t2 are too close to each other.");
  else if (status == CV_TOO_MUCH_WORK)
  {
    char err[1024];
    snprintf(err, 1024, "solver stopped at t = %g after maximum number of steps.", current_time);
    status_message = string_dup(err);
  }
  else if (status == CV_TOO_MUCH_ACC)
    status_message = string_dup("solver could not achieve desired level of accuracy.");
  else if (status == CV_ERR_FAILURE)
    status_message = string_dup("solver encountered too many error test failures.");
  else if (status == CV_CONV_FAILURE)
    status_message = string_dup("solver encountered too many convergence test failures.");
  else if (status == CV_LINIT_FAIL)
    status_message = string_dup("solver's linear solver failed to initialize.");
  else if (status == CV_LSETUP_FAIL)
    status_message = string_dup("solver's linear solver setup failed.");
  else if (status == CV_LSOLVE_FAIL)
    status_message = string_dup("solver's linear solver failed.");
  else if (status == CV_RHSFUNC_FAIL)
    status_message = string_dup("solver's RHS function failed unrecoverably.");
  else if (status == CV_FIRST_RHSFUNC_ERR)
    status_message = string_dup("solver's first call to RHS function failed.");
  else if (status == CV_REPTD_RHSFUNC_ERR)
    status_message = string_dup("solver encountered too many recoverable RHS failures.");
  else if (status == CV_UNREC_RHSFUNC_ERR)
    status_message = string_dup("solver failed to recover from a recoverable RHS failure.");
  else if (status == CV_RTFUNC_FAIL)
    status_message = string_dup("solver encountered a failure in the rootfinding function.");
  return status_message;
}

// This function wraps around the user-supplied right hand side.
static int bdf_evaluate_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  START_FUNCTION_TIMER();
  bdf_ode_t* solver = context;
  real_t* xx = NV_DATA(U);
  real_t* xxd = NV_DATA(U_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(solver->U_with_ghosts, xx, sizeof(real_t) * solver->num_local_values);
  int status = solver->rhs(solver->context, t, solver->U_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, solver->U_with_ghosts, sizeof(real_t) * solver->num_local_values);
  }

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < solver->observers->size; ++i)
  {
    bdf_ode_observer_t* obs = solver->observers->data[i];
    if (obs->rhs_computed != NULL)
      obs->rhs_computed(obs->context, t, xx, xxd);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

static bool bdf_step(void* context, real_t max_dt, real_t* t, real_t* U)
{
  START_FUNCTION_TIMER();
  bdf_ode_t* solver = context;
  int status = CV_SUCCESS;

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > solver->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = CVode(solver->cvode, t2, solver->U, &solver->t, CV_ONE_STEP);
    if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN))
    {
      solver->status_message = get_status_message(status, solver->t);
      STOP_FUNCTION_TIMER();
      return false;
    }

    if ((t2 - *t) < (solver->t - *t))
      log_detail("bdf_ode_solver: took internal step dt = %g", solver->t - *t);
  }

  // If we integrated past t2, interpolate to t2.
  if (solver->t > t2)
  {
    status = CVodeGetDky(solver->cvode, t2, 0, solver->U);
    *t = t2;
  }
  else
    *t = solver->t;
  
  // Clear the present status.
  if (solver->status_message != NULL)
  {
    polymec_free(solver->status_message);
    solver->status_message = NULL;
  }

  // Did it work?
  if ((status == CV_SUCCESS) || (status == CV_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(solver->U), sizeof(real_t) * solver->num_local_values); 
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    solver->status_message = get_status_message(status, solver->t);
    STOP_FUNCTION_TIMER();
    return false;
  }
}

static bool bdf_advance(void* context, real_t t1, real_t t2, real_t* U)
{
  bdf_ode_t* solver = context;
  solver->t = t1;
  CVodeReInit(solver->cvode, t1, solver->U);
  CVodeSetStopTime(solver->cvode, t2);
  
  // Copy in the solution.
  memcpy(NV_DATA(solver->U), U, sizeof(real_t) * solver->num_local_values); 

  // Integrate.
  int status = CVode(solver->cvode, t2, solver->U, &solver->t, CV_NORMAL);
  
  // Clear the present status.
  if (solver->status_message != NULL)
  {
    polymec_free(solver->status_message);
    solver->status_message = NULL;
  }

  // Did it work?
  if ((status == CV_SUCCESS) || (status == CV_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(solver->U), sizeof(real_t) * solver->num_local_values); 
    return true;
  }
  else
  {
    solver->status_message = get_status_message(status, solver->t);
    return false;
  }
}

static void bdf_reset(void* context, real_t t, real_t* U)
{
  bdf_ode_t* solver = context;

  // Reset the preconditioner.
  if (solver->precond != NULL)
    newton_pc_reset(solver->precond, t);

  // Copy in the solution and reinitialize.
  memcpy(NV_DATA(solver->U), U, sizeof(real_t) * solver->num_local_values); 
  CVodeReInit(solver->cvode, t, solver->U);
  solver->t = t;
}

static void bdf_dtor(void* context)
{
  bdf_ode_t* solver = context;

  // Kill the preconditioner stuff.
  if (solver->precond != NULL)
    newton_pc_free(solver->precond);

  // Kill the CVode stuff.
  polymec_free(solver->U_with_ghosts);
  N_VDestroy(solver->U);
  if (solver->nls != NULL)
    SUNNonlinSolFree(solver->nls);
  if (solver->ls != NULL)
    SUNLinSolFree(solver->ls);
  CVodeFree(&solver->cvode);

  // Kill the rest.
  if (solver->status_message != NULL)
    polymec_free(solver->status_message);
  if ((solver->context != NULL) && (solver->dtor != NULL))
    solver->dtor(solver->context);
  ptr_array_free(solver->observers);
  if (solver->error_weights != NULL)
    polymec_free(solver->error_weights);
  polymec_free(solver);
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(real_t t, N_Vector U, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context)
{
  bdf_ode_t* solver = context;
  if (!jacobian_is_current)
  {
    // Compute the approximate Jacobian using a solution vector with ghosts.
    log_debug("jfnk_bdf_ode_solver: Calculating P = I - %g * J.", gamma);
    memcpy(solver->U_with_ghosts, NV_DATA(U), sizeof(real_t) * solver->num_local_values);
    newton_pc_setup(solver->precond, 1.0, -gamma, 0.0, t, solver->U_with_ghosts, NULL);
    *jacobian_was_updated = 1;
  }
  else
    *jacobian_was_updated = 0;
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector U, N_Vector F, 
                                       N_Vector r, N_Vector z, 
                                       real_t gamma, real_t delta, 
                                       int lr, void* context)
{
  bdf_ode_t* solver = context;
  
  // FIXME: Apply scaling if needed.

  // Set the preconditioner's L2 norm tolerance.
  newton_pc_set_tolerance(solver->precond, delta);

  // Solve it.
  int result = (newton_pc_solve(solver->precond, t, NV_DATA(U), NULL,
                                NV_DATA(r), NV_DATA(z))) ? 0 : 1;
  return result;
}

// Adaptor for J*y setup function
static int set_up_Jy(real_t t, N_Vector y, N_Vector Jy, void* context)
{
  // Nothing here!
  return 0;
}

// Adaptor for J*y function
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector U, N_Vector rhs, void* context, N_Vector tmp)
{
  START_FUNCTION_TIMER();
  bdf_ode_t* solver = context;
  real_t* xx = NV_DATA(U);
  real_t* yy = NV_DATA(y);
  real_t* my_rhs = NV_DATA(rhs);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);

  // Make sure we use ghosts.
  memcpy(solver->U_with_ghosts, xx, sizeof(real_t) * solver->num_local_values);
  int status = solver->Jy(solver->context, t, xx, my_rhs, yy, temp, jjy);

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < solver->observers->size; ++i)
  {
    bdf_ode_observer_t* obs = solver->observers->data[i];
    if (obs->Jy_computed != NULL)
      obs->Jy_computed(obs->context, t, xx, my_rhs, yy, jjy);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

ode_solver_t* jfnk_bdf_ode_solver_new(int order,
                                      MPI_Comm comm,
                                      int num_local_values, 
                                      int num_remote_values, 
                                      void* context, 
                                      int (*rhs_func)(void* context, real_t t, real_t* U, real_t* xdot),
                                      int (*Jy_func)(void* context, real_t t, real_t* U, real_t* rhs, real_t* y, real_t* Jy, real_t* temp),
                                      void (*dtor)(void* context),
                                      newton_pc_t* precond,
                                      jfnk_bdf_krylov_t solver_type,
                                      int max_krylov_dim)
{
  ASSERT(order >= 1);
  ASSERT(order <= 5);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs_func != NULL);
  ASSERT(precond != NULL);
  ASSERT(max_krylov_dim > 3);
  ASSERT(!newton_pc_coefficients_fixed(precond));

  bdf_ode_t* solver = polymec_malloc(sizeof(bdf_ode_t));
  solver->comm = comm;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->context = context;
  solver->rhs = rhs_func;
  solver->dtor = dtor;
  solver->status_message = NULL;
  solver->max_krylov_dim = max_krylov_dim;
  solver->Jy = Jy_func;
  solver->t = 0.0;
  solver->observers = ptr_array_new();
  solver->error_weights = NULL;

  solver->reset_func = NULL;
  solver->setup_func = NULL;
  solver->solve_func = NULL;

  // Set up CVode and accessories.
  solver->U = N_VNew(solver->comm, solver->num_local_values);
  solver->U_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->cvode = CVodeCreate(CV_BDF);
  CVodeSetMaxOrd(solver->cvode, order);
  CVodeSetUserData(solver->cvode, solver);
  CVodeInit(solver->cvode, bdf_evaluate_rhs, 0.0, solver->U);
  solver->nls = SUNNonlinSol_Newton(solver->U);
  CVodeSetNonlinearSolver(solver->cvode, solver->nls);

  newton_pc_side_t newton_side = newton_pc_side(precond);
  int side;
  switch (newton_side)
  {
    case NEWTON_PC_LEFT:
      side = PREC_LEFT;
      break;
    case NEWTON_PC_RIGHT:
      side = PREC_RIGHT;
      break;
    case NEWTON_PC_BOTH:
      side = PREC_BOTH;
  }

  // Set up the solver type.
  if (solver_type == JFNK_BDF_GMRES)
  {
    solver->ls = SUNSPGMR(solver->U, side, max_krylov_dim);
    // We use modified Gram-Schmidt orthogonalization.
    SUNSPGMRSetGSType(solver->ls, MODIFIED_GS);
    CVodeSetLinearSolver(solver->cvode, solver->ls, NULL);
  }
  else if (solver_type == JFNK_BDF_BICGSTAB)
  {
    solver->ls = SUNSPBCGS(solver->U, side, max_krylov_dim);
    CVodeSetLinearSolver(solver->cvode, solver->ls, NULL);
  }
  else
  {
    solver->ls = SUNSPTFQMR(solver->U, side, max_krylov_dim);
    CVodeSetLinearSolver(solver->cvode, solver->ls, NULL);
  }

  // Set up the Jacobian function and preconditioner.
  if (Jy_func != NULL)
    CVodeSetJacTimes(solver->cvode, set_up_Jy, eval_Jy);
  solver->precond = precond;
  CVodeSetPreconditioner(solver->cvode, set_up_preconditioner,
                           solve_preconditioner_system);

  ode_solver_vtable vtable = {.step = bdf_step, .advance = bdf_advance, .reset = bdf_reset, .dtor = bdf_dtor};
  char name[1024];
  snprintf(name, 1024, "JFNK Backwards-Difference-Formulae (order %d)", order);
  ode_solver_t* I = ode_solver_new(name, solver, vtable, order,
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_solver_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void* bdf_ode_solver_context(ode_solver_t* solver)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  return bdf->context;
}

void bdf_ode_solver_set_max_err_test_failures(ode_solver_t* solver,
                                              int max_failures)
{
  ASSERT(max_failures > 0);
  bdf_ode_t* bdf = ode_solver_context(solver);
  CVodeSetMaxErrTestFails(bdf->cvode, max_failures);
}

void bdf_ode_solver_set_max_nonlinear_iterations(ode_solver_t* solver,
                                                 int max_iterations)
{
  ASSERT(max_iterations > 0);
  bdf_ode_t* bdf = ode_solver_context(solver);
  CVodeSetMaxNonlinIters(bdf->cvode, max_iterations);
}

void bdf_ode_solver_set_nonlinear_convergence_coeff(ode_solver_t* solver,
                                                    real_t coefficient)
{
  ASSERT(coefficient > 0.0);
  bdf_ode_t* bdf = ode_solver_context(solver);
  CVodeSetNonlinConvCoef(bdf->cvode, (double)coefficient);
}

void bdf_ode_solver_set_tolerances(ode_solver_t* solver,
                                      real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  bdf_ode_t* bdf = ode_solver_context(solver);

  // Clear any existing error weight function.
  bdf->compute_weights = NULL;
  if (bdf->error_weights != NULL)
  {
    polymec_free(bdf->error_weights);
    bdf->error_weights = NULL;
  }

  // Set the tolerances.
  CVodeSStolerances(bdf->cvode, relative_tol, absolute_tol);
}

// Constant error weight adaptor function.
static void use_constant_weights(void* context, real_t* y, real_t* weights)
{
  bdf_ode_t* solver = context;
  ASSERT(solver->error_weights != NULL);
  memcpy(weights, solver->error_weights, sizeof(real_t) * solver->num_local_values);
}

void bdf_ode_solver_set_error_weights(ode_solver_t* solver, real_t* weights)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
#ifndef NDEBUG
  // Check for non-negativity and total positivity.
  real_t total = 0.0;
  for (int i = 0; i < bdf->num_local_values; ++i)
  {
    ASSERT(weights[i] >= 0.0);
    total += weights[i];
  }
  ASSERT(total > 0.0);
#endif

  if (bdf->error_weights == NULL)
    bdf->error_weights = polymec_malloc(sizeof(real_t) * bdf->num_local_values);
  memcpy(bdf->error_weights, weights, sizeof(real_t) * bdf->num_local_values);
  bdf_ode_solver_set_error_weight_function(solver, use_constant_weights);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  bdf_ode_t* solver = context;
  solver->compute_weights(solver->context, NV_DATA(y), NV_DATA(ewt));

  // Check that all the weights are non-negative.
  int N = (int)NV_LOCLENGTH(y);
  for (int i = 0; i < N; ++i)
  {
    if (NV_Ith(y, i) < 0.0)
      return -1;
  }
  return 0;
}

void bdf_ode_solver_set_error_weight_function(ode_solver_t* solver,
                                              void (*compute_weights)(void* context, real_t* y, real_t* weights))
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  ASSERT(compute_weights != NULL);
  bdf->compute_weights = compute_weights;
  CVodeWFtolerances(bdf->cvode, compute_error_weights);
}

void bdf_ode_solver_set_stability_limit_detection(ode_solver_t* solver,
                                                  bool use_detection)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  CVodeSetStabLimDet(bdf->cvode, use_detection);
}

void bdf_ode_solver_eval_rhs(ode_solver_t* solver, real_t t, real_t* U, real_t* U_dot)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  memcpy(bdf->U_with_ghosts, U, sizeof(real_t) * bdf->num_local_values);
  bdf->rhs(bdf->context, t, bdf->U_with_ghosts, U_dot);
}

newton_pc_t* bdf_ode_solver_preconditioner(ode_solver_t* solver)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  return bdf->precond;
}

void bdf_ode_solver_get_diagnostics(ode_solver_t* solver, 
                                    bdf_ode_solver_diagnostics_t* diagnostics)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  diagnostics->status_message = bdf->status_message; // borrowed!
  CVodeGetNumSteps(bdf->cvode, &diagnostics->num_steps);
  CVodeGetLastOrder(bdf->cvode, &diagnostics->order_of_last_step);
  CVodeGetCurrentOrder(bdf->cvode, &diagnostics->order_of_next_step);
  CVodeGetLastStep(bdf->cvode, &diagnostics->last_step_size);
  CVodeGetCurrentStep(bdf->cvode, &diagnostics->next_step_size);
  CVodeGetNumRhsEvals(bdf->cvode, &diagnostics->num_rhs_evaluations);
  CVodeGetNumLinSolvSetups(bdf->cvode, &diagnostics->num_linear_solve_setups);
  CVodeGetNumErrTestFails(bdf->cvode, &diagnostics->num_error_test_failures);
  CVodeGetNumNonlinSolvIters(bdf->cvode, &diagnostics->num_nonlinear_solve_iterations);
  CVodeGetNumNonlinSolvConvFails(bdf->cvode, &diagnostics->num_nonlinear_solve_convergence_failures);
  if (bdf->solve_func == NULL) // JFNK mode
  {
    CVodeGetNumLinIters(bdf->cvode, &diagnostics->num_linear_solve_iterations);
    CVodeGetNumPrecEvals(bdf->cvode, &diagnostics->num_preconditioner_evaluations);
    CVodeGetNumPrecSolves(bdf->cvode, &diagnostics->num_preconditioner_solves);
    CVodeGetNumLinConvFails(bdf->cvode, &diagnostics->num_linear_solve_convergence_failures);
  }
  else
  {
    diagnostics->num_linear_solve_iterations = (long int)bdf->num_linear_iterations;
    diagnostics->num_linear_solve_convergence_failures = (long int)bdf->num_linear_conv_failures;
    diagnostics->num_preconditioner_solves = -1;
    diagnostics->num_preconditioner_evaluations = -1;
  }
}

void bdf_ode_solver_diagnostics_fprintf(bdf_ode_solver_diagnostics_t* diagnostics, 
                                        FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "ODE solver diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Order of last step: %d\n", diagnostics->order_of_last_step);
  fprintf(stream, "  Order of next step: %d\n", diagnostics->order_of_next_step);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Next step size: %g\n", diagnostics->next_step_size);
  fprintf(stream, "  Num RHS evaluations: %d\n", (int)diagnostics->num_rhs_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  if (diagnostics->num_preconditioner_evaluations != -1) // JFNK mode
  {
    fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
    fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
  }
}

bdf_ode_observer_t* bdf_ode_observer_new(void* context,
                                         void (*rhs_computed)(void* context, real_t t, real_t* U, real_t* rhs),
                                         void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* rhs, real_t* y, real_t* Jy),
                                         void (*dtor)(void* context))
{
  bdf_ode_observer_t* obs = polymec_malloc(sizeof(bdf_ode_observer_t));
  obs->context = context;
  obs->rhs_computed = rhs_computed;
  obs->Jy_computed = Jy_computed;
  obs->dtor = dtor;
  return obs;
}

static void bdf_ode_observer_free(bdf_ode_observer_t* observer)
{
  if ((observer->dtor != NULL) && (observer->context != NULL))
    observer->dtor(observer->context);
  polymec_free(observer);
}

void bdf_ode_solver_add_observer(ode_solver_t* solver,
                                 bdf_ode_observer_t* observer)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  ptr_array_append_with_dtor(bdf->observers, observer, DTOR(bdf_ode_observer_free));
}

// This unpublished function provides direct access to the CVode object 
// within the BDF solver.
void* bdf_ode_solver_cvode(ode_solver_t* solver);
void* bdf_ode_solver_cvode(ode_solver_t* solver)
{
  bdf_ode_t* bdf = ode_solver_context(solver);
  return bdf->cvode;
}

// This unpublished function provides direct access to the CVode right-hand
// side function within the BDF solver.
int bdf_ode_solver_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context);
int bdf_ode_solver_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  return bdf_evaluate_rhs(t, U, U_dot, context);
}

//------------------------------------------------------------------------
//                    Custom BDF solver stuff
//------------------------------------------------------------------------

static int bdf_linit(CVodeMem cv_mem)
{
  bdf_ode_t* bdf = cv_mem->cv_user_data;
  real_t t = cv_mem->cv_tn;
  real_t* U = NV_DATA(cv_mem->cv_y);
  bdf->sqrtN = -1.0;
  bdf->num_linear_iterations = 0;
  bdf->num_linear_conv_failures = 0;
  return bdf->reset_func(bdf->context, t, U);
}

static int bdf_lsetup(CVodeMem cv_mem, 
                      int convfail,
                      N_Vector ypred,
                      N_Vector fpred,
                      booleantype* jcurPtr,
                      N_Vector vtemp1,
                      N_Vector vtemp2,
                      N_Vector vtemp3)
{
  bdf_ode_t* bdf = cv_mem->cv_user_data;
  bdf_conv_status_t conv_status;
  if (convfail == CV_NO_FAILURES)
    conv_status = BDF_CONV_NO_FAILURES;
  else if (convfail == CV_FAIL_BAD_J)
    conv_status = BDF_CONV_BAD_J_FAILURE;
  else
    conv_status = BDF_CONV_OTHER_FAILURE;
  real_t gamma = cv_mem->cv_gamma;
  int step = (int)cv_mem->cv_nst;
  real_t t = cv_mem->cv_tn;
  bool J_updated = false;
  real_t* U_pred = NV_DATA(ypred);
  real_t* U_dot_pred = NV_DATA(fpred);
  real_t* work1 = NV_DATA(vtemp1);
  real_t* work2 = NV_DATA(vtemp2);
  real_t* work3 = NV_DATA(vtemp3);
  if (bdf->sqrtN <= 0.0)
  {
    N_VConst(1.0, vtemp1);
    bdf->sqrtN = sqrt(N_VDotProd(vtemp1, vtemp1));
  }
  int status = bdf->setup_func(bdf->context, conv_status, gamma, step, t, 
                               U_pred, U_dot_pred, &J_updated, work1, 
                               work2, work3);
  *jcurPtr = J_updated;
  return status;
}

static int bdf_lsolve(CVodeMem cv_mem, 
                      N_Vector b, 
                      N_Vector weight,
                      N_Vector ycur, 
                      N_Vector fcur)
{
  bdf_ode_t* bdf = cv_mem->cv_user_data;
  real_t t = cv_mem->cv_tn;
  real_t* U = NV_DATA(ycur);
  real_t* U_dot = NV_DATA(fcur);
  real_t* W = NV_DATA(weight);
  real_t C1_inv = cv_mem->cv_tq[4]; // constant 1/C' in nonlinear iteration conv test.
//  real_t res_norm_tol = 0.05 * bdf->sqrtN * C1_inv; // FIXME: Why doesn't this work?
  real_t res_norm_tol = 0.05 * C1_inv;
  real_t* B = NV_DATA(b);
  int num_iters = 0;
  int result = bdf->solve_func(bdf->context, t, U, U_dot, W, res_norm_tol, B, &num_iters);
  bdf->num_linear_iterations += num_iters;
  if (result != 0)
    ++(bdf->num_linear_conv_failures);
  return result;
}

static int bdf_lfree(struct CVodeMemRec* cv_mem)
{
  return 0;
}

ode_solver_t* bdf_ode_solver_new(const char* name,
                                 int order, 
                                 MPI_Comm comm,
                                 int num_local_values, 
                                 int num_remote_values, 
                                 void* context, 
                                 int (*rhs_func)(void* context, real_t t, real_t* U, real_t* xdot),
                                 int (*reset_func)(void* context, real_t t, real_t* U),
                                 int (*setup_func)(void* context, 
                                                   bdf_conv_status_t conv_status, 
                                                   real_t gamma, 
                                                   int step,
                                                   real_t t, 
                                                   real_t* U_pred, 
                                                   real_t* U_dot_pred, 
                                                   bool* J_current, 
                                                   real_t* work1, real_t* work2, real_t* work3),
                                 int (*solve_func)(void* context, 
                                                   real_t t, 
                                                   real_t* U,
                                                   real_t* U_dot,
                                                   real_t* W, 
                                                   real_t res_norm_tol, 
                                                   real_t* B,
                                                   int* num_iters),
                                 void (*dtor)(void* context))
{
  ASSERT(order >= 1);
  ASSERT(order <= 5);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs_func != NULL);
  ASSERT(reset_func != NULL);
  ASSERT(setup_func != NULL);
  ASSERT(solve_func != NULL);

  bdf_ode_t* solver = polymec_malloc(sizeof(bdf_ode_t));
  solver->comm = comm;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->context = context;
  solver->rhs = rhs_func;
  solver->dtor = dtor;
  solver->status_message = NULL;
  solver->max_krylov_dim = -1;
  solver->Jy = NULL;
  solver->precond = NULL;
  solver->t = 0.0;
  solver->observers = ptr_array_new();
  solver->error_weights = NULL;

  // Set up CVode and accessories.
  solver->U = N_VNew(solver->comm, solver->num_local_values);
  solver->U_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->cvode = CVodeCreate(CV_BDF);
  CVodeSetMaxOrd(solver->cvode, order);
  CVodeSetUserData(solver->cvode, solver);
  CVodeInit(solver->cvode, bdf_evaluate_rhs, 0.0, solver->U);
  solver->nls = SUNNonlinSol_Newton(solver->U);
  CVodeSetNonlinearSolver(solver->cvode, solver->nls);

  // Set up the solver.
  solver->ls = NULL;
  solver->reset_func = reset_func;
  solver->setup_func = setup_func;
  solver->solve_func = solve_func;
  CVodeMem cv_mem = solver->cvode;
  cv_mem->cv_linit = bdf_linit;
  cv_mem->cv_lsetup = bdf_lsetup;
  cv_mem->cv_lsolve = bdf_lsolve;
  cv_mem->cv_lfree = bdf_lfree;

  ode_solver_vtable vtable = {.step = bdf_step, 
                                  .advance = bdf_advance, 
                                  .reset = bdf_reset, 
                                  .dtor = bdf_dtor};
  ode_solver_t* I = ode_solver_new(name, solver, vtable, order,
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_solver_set_tolerances(I, 1e-4, 1.0);

  return I;
}

//------------------------------------------------------------------------
//                  Inexact Newton-Krylov solver stuff
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;

  // Behavior.
  void* context;
  int (*rhs_func)(void* context, real_t t, real_t* U, real_t* U_dot);
  int (*J_func)(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J);
  void (*dtor)(void* context);

  // Linear system stuff.
  krylov_factory_t* factory;
  krylov_solver_t* solver;
  krylov_pc_t* pc;
  krylov_matrix_t* M; // Newton matrix M = (I - gamma * J).
  krylov_vector_t* X; // Solution vector.
  krylov_vector_t* B; // Right-hand side vector.
  krylov_vector_t* W; // Weights vector.

  // Metadata.
  real_t gamma_prev; // Previous value of gamma in Newton matrix.
  int prev_J_update_step; // Step during which J was previously updated.
  matrix_sparsity_t* sparsity;
  int block_size;
} ink_bdf_ode_t;

static int ink_rhs(void* context, real_t t, real_t* U, real_t* U_dot)
{
  ink_bdf_ode_t* ink = context;
  return ink->rhs_func(ink->context, t, U, U_dot);
}

static int ink_reset(void* context, real_t t, real_t* U)
{
  START_FUNCTION_TIMER();
  ink_bdf_ode_t* ink = context;
  
  // Free any resources currently in use.
  if (ink->M != NULL)
    krylov_matrix_free(ink->M);
  if (ink->block_size > 1)
    ink->M = krylov_factory_block_matrix(ink->factory, ink->sparsity, ink->block_size);
  else
    ink->M = krylov_factory_matrix(ink->factory, ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);

  // Allocate resources.
  index_t* row_dist = matrix_sparsity_row_distribution(ink->sparsity);
  ink->X = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->B = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->W = krylov_factory_vector(ink->factory, ink->comm, row_dist);

  // Set some metadata.
  ink->gamma_prev = 1.0;
  ink->prev_J_update_step = 0;

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_setup(void* context, 
                     bdf_conv_status_t conv_status, 
                     real_t gamma, 
                     int step,
                     real_t t, 
                     real_t* U_pred, 
                     real_t* U_dot_pred, 
                     bool* J_updated, 
                     real_t* work1, real_t* work2, real_t* work3)
{
  START_FUNCTION_TIMER();
  ink_bdf_ode_t* ink = context;

  real_t dgamma = ABS(gamma / ink->gamma_prev - 1.0);
  log_debug("ink_bdf_ode_solver: Calculating M = I - %g * J.", gamma);
  if ((step == 0) ||
      (step > ink->prev_J_update_step + 50) ||
      ((conv_status == BDF_CONV_BAD_J_FAILURE) && (dgamma < 0.2)) || 
      (conv_status == BDF_CONV_OTHER_FAILURE))
  {
    char reason[129];
    // Call our Jacobian calculation function.
    if (step == 0)
      snprintf(reason, 128, "first step");
    else if (step > ink->prev_J_update_step + 50)
      snprintf(reason, 128, "> 50 steps since last calculation");
    else if (conv_status == BDF_CONV_BAD_J_FAILURE)
      snprintf(reason, 128, "outdated Newton matrix");
    else
      snprintf(reason, 128, "convergence failure reduced dt");
    log_debug("ink_bdf_ode_solver: Updating J (reason: %s).", reason);
    int status = ink->J_func(ink->context, t, U_pred, U_dot_pred, ink->M);
    if (status != 0)
      return status;

    // Scale the Jacobian by -gamma and add the identity matrix.
    krylov_matrix_scale(ink->M, -gamma);
    krylov_matrix_add_identity(ink->M, 1.0);

    // Use this matrix as the operator in our solver.
    krylov_solver_set_operator(ink->solver, ink->M);

    // Save some information.
    ink->prev_J_update_step = step;
    ink->gamma_prev = gamma;
    *J_updated = true;
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    if (conv_status == BDF_CONV_NO_FAILURES)
    {
      // We don't need to update J, but we still need to recompute I - gamma*J.
      krylov_matrix_add_identity(ink->M, -1.0);
      krylov_matrix_scale(ink->M, gamma/ink->gamma_prev);
      krylov_matrix_add_identity(ink->M, 1.0);
      ink->gamma_prev = gamma;
      *J_updated = false;
      STOP_FUNCTION_TIMER();
      return 0;
    }
    else
    {
      STOP_FUNCTION_TIMER();
      return 1;
    }
  }
}

static int ink_solve(void* context, 
                     real_t t, 
                     real_t* U,
                     real_t* U_dot,
                     real_t* W, 
                     real_t res_norm_tol,
                     real_t* B, 
                     int* num_iters) 
{
  START_FUNCTION_TIMER();
  ink_bdf_ode_t* ink = context;

  // Copy RHS data from B into ink->B.
  krylov_vector_copy_in(ink->B, B);

  // Copy weights into ink->W.
  krylov_vector_copy_in(ink->W, W);

  // If the WRMS norm of B is less than our tolerance, return X = 0.
  real_t B_norm = krylov_vector_wrms_norm(ink->B, ink->W);
  if (B_norm < res_norm_tol)
  {
    log_debug("ink_bdf_ode_solver: ||B|| < tolerance (%g < %g), so X -> 0.", B_norm, res_norm_tol); 
    krylov_vector_zero(ink->X);
    krylov_vector_copy_out(ink->X, B);
    STOP_FUNCTION_TIMER();
    return 0;
  }

  // Set the tolerance on the residual norm.
  real_t rel_tol = 1e-8;
  real_t div_tol = 10.0;
  krylov_solver_set_tolerances(ink->solver, rel_tol, res_norm_tol, div_tol);

  // Solve A*X = B.
  real_t res_norm;
  bool solved = krylov_solver_solve_scaled(ink->solver, ink->B, ink->W, ink->W, 
                                           ink->X, &res_norm, num_iters);

  if (solved)
  {
    log_debug("ink_bdf_ode_solver: Solved A*X = B (||R|| == %g after %d iters).", res_norm, num_iters);

    // Copy solution data from ink->X into B.
    krylov_vector_copy_out(ink->X, B);
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    log_debug("ink_bdf_ode_solver: Solution to A*X = B did not converge.");
    STOP_FUNCTION_TIMER();
    return 1;
  }
}

static void ink_dtor(void* context)
{
  ink_bdf_ode_t* ink = context;
  matrix_sparsity_free(ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);
  if (ink->W != NULL)
    krylov_vector_free(ink->W);
  if (ink->M != NULL)
    krylov_matrix_free(ink->M);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  krylov_solver_free(ink->solver);
  krylov_factory_free(ink->factory);
  if ((ink->context != NULL) && (ink->dtor != NULL))
    ink->dtor(ink->context);
}

ode_solver_t* ink_bdf_ode_solver_new(int order, 
                                     MPI_Comm comm,
                                     krylov_factory_t* factory,
                                     matrix_sparsity_t* J_sparsity,
                                     void* context, 
                                     int (*rhs_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                                     int (*J_func)(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J),
                                     void (*dtor)(void* context))
{
  ink_bdf_ode_t* ink = polymec_malloc(sizeof(ink_bdf_ode_t));
  ink->comm = comm;
  ink->context = context;
  ink->rhs_func = rhs_func;
  ink->J_func = J_func;
  ink->dtor = dtor;
  ink->factory = factory;
  ink->sparsity = J_sparsity;
  ink->solver = krylov_factory_gmres_solver(ink->factory, comm, 30);
  ink->pc = NULL;
  ink->M = NULL;
  ink->B = NULL;
  ink->X = NULL;
  ink->W = NULL;
  ink->block_size = 1;

  char name[1024];
  snprintf(name, 1024, "INK Backwards-Difference-Formulae (order %d)", order);
  int num_local_values = (int)(matrix_sparsity_num_local_rows(J_sparsity));
  ode_solver_t* I = bdf_ode_solver_new(name, order, comm, 
                                       num_local_values, 0,
                                       ink, ink_rhs, ink_reset, 
                                       ink_setup, ink_solve, ink_dtor);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_solver_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void ink_bdf_ode_solver_use_pcg(ode_solver_t* ink_bdf_ode_solver)
{
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_pcg_solver(ink->factory, ink->comm);
}

void ink_bdf_ode_solver_use_gmres(ode_solver_t* ink_bdf_ode_solver,
                                  int max_krylov_dim)
{
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_gmres_solver(ink->factory, ink->comm, max_krylov_dim);
}

void ink_bdf_ode_solver_use_bicgstab(ode_solver_t* ink_bdf_ode_solver)
{
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_bicgstab_solver(ink->factory, ink->comm);
}

void ink_bdf_ode_solver_use_special(ode_solver_t* ink_bdf_ode_solver,
                                    const char* solver_name,
                                    string_string_unordered_map_t* options)
{
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_special_solver(ink->factory, ink->comm,
                                              solver_name, options);
}

void ink_bdf_ode_solver_set_pc(ode_solver_t* ink_bdf_ode_solver,
                               const char* pc_name, 
                               string_string_unordered_map_t* options)
{
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  ink->pc = krylov_factory_preconditioner(ink->factory, ink->comm, pc_name, options);
}

void ink_bdf_ode_solver_set_block_size(ode_solver_t* ink_bdf_ode_solver,
                                       int block_size)
{
  ASSERT(block_size > 0);
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  ink->block_size = block_size;
}

void* ink_bdf_ode_solver_context(ode_solver_t* ink_bdf_ode_solver)
{
  ink_bdf_ode_t* ink = bdf_ode_solver_context(ink_bdf_ode_solver);
  return ink->context;
}
