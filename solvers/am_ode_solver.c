// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "core/array.h"
#include "core/timer.h"
#include "solvers/am_ode_solver.h"
#include "solvers/newton_pc.h"
#include "cvode/cvode.h"

// JFNK stuff.
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunlinsol/sunlinsol_spfgmr.h"
#include "sunlinsol/sunlinsol_spbcgs.h"
#include "sunlinsol/sunlinsol_sptfqmr.h"

#include "sunnonlinsol/sunnonlinsol_fixedpoint.h"
#include "sunnonlinsol/sunnonlinsol_newton.h"

struct am_ode_observer_t 
{
  void* context;
  void (*rhs_computed)(void* context, real_t t, real_t* U, real_t* rhs);
  void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* rhs, real_t* y, real_t* Jy);
  void (*dtor)(void* context);
};

typedef struct
{
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context; 
  int (*rhs)(void* context, real_t t, real_t* U, real_t* U_dot);
  void (*dtor)(void* context);

  // CVODE data structures.
  void* cvode;
  N_Vector U; 
  real_t* U_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  int (*Jy)(void* context, real_t t, real_t* U, real_t* rhs, real_t* y, real_t* temp, real_t* Jy);

  // Error weight function.
  void (*compute_weights)(void* context, real_t* y, real_t* weights);
  real_t* error_weights;

  // Observers.
  ptr_array_t* observers;
} am_ode_t;

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
static int am_evaluate_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  START_FUNCTION_TIMER();
  am_ode_t* solver = context;
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
    am_ode_observer_t* obs = solver->observers->data[i];
    if (obs->rhs_computed != NULL)
      obs->rhs_computed(obs->context, t, xx, xxd);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

static bool am_step(void* context, real_t max_dt, real_t* t, real_t* U)
{
  START_FUNCTION_TIMER();
  am_ode_t* solver = context;
  int status = CV_SUCCESS;

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > solver->t)
  {
    // solver to at least t -> t + max_dt.
    status = CVode(solver->cvode, t2, solver->U, &solver->t, CV_ONE_STEP);
    if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN))
    {
      solver->status_message = get_status_message(status, solver->t);
      STOP_FUNCTION_TIMER();
      return false;
    }
    if ((t2 - *t) < (solver->t - *t))
      log_detail("am_ode_solver: took internal step dt = %g", solver->t - *t);
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

static bool am_advance(void* context, real_t t1, real_t t2, real_t* U)
{
  am_ode_t* solver = context;
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

static void am_reset(void* context, real_t t, real_t* U)
{
  am_ode_t* solver = context;

  // Reset the preconditioner.
  if (solver->precond != NULL)
    newton_pc_reset(solver->precond, t);

  // Copy in the solution and reinitialize.
  memcpy(NV_DATA(solver->U), U, sizeof(real_t) * solver->num_local_values); 
  CVodeReInit(solver->cvode, t, solver->U);
  solver->t = t;
}

static void am_dtor(void* context)
{
  am_ode_t* solver = context;

  // Kill the preconditioner stuff.
  if (solver->precond != NULL)
    newton_pc_free(solver->precond);

  // Kill the CVode stuff.
  polymec_free(solver->U_with_ghosts);
  N_VDestroy(solver->U);
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

ode_solver_t* functional_am_ode_solver_new(int order, 
                                           MPI_Comm comm, 
                                           int num_local_values,
                                           int num_remote_values,
                                           void* context, 
                                           int (*rhs)(void* context, real_t t, real_t* U, real_t* U_dot),
                                           void (*dtor)(void* context))
{
  ASSERT(order >= 1);
  ASSERT(order <= 12);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs != NULL);

  am_ode_t* solver = polymec_malloc(sizeof(am_ode_t));
  solver->comm = comm;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->context = context;
  solver->rhs = rhs;
  solver->dtor = dtor;
  solver->status_message = NULL;
  solver->max_krylov_dim = 0;
  solver->Jy = NULL;
  solver->precond = NULL;
  solver->observers = ptr_array_new();
  solver->error_weights = NULL;

  // Set up KINSol and accessories.
  solver->U = N_VNew(solver->comm, solver->num_local_values);
  solver->U_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->cvode = CVodeCreate(CV_ADAMS);
  CVodeSetMaxOrd(solver->cvode, order);
  CVodeSetUserData(solver->cvode, solver);
  CVodeInit(solver->cvode, am_evaluate_rhs, 0.0, solver->U);
  int num_accels; // FIXME: Number of Anderson acceleration vectors!
  SUNNonlinearSolver nls = SUNNonlinSol_FixedPoint(solver->U, num_accels);
  CVodeSetNonlinearSolver(solver->cvode, nls);

  ode_solver_vtable vtable = {.step = am_step, .advance = am_advance, .reset = am_reset, .dtor = am_dtor};
  char name[1024];
  snprintf(name, 1024, "Functional Adams-Moulton (order %d)", order);
  ode_solver_t* I = ode_solver_new(name, solver, vtable, order,
                                   num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  am_ode_solver_set_tolerances(I, 1e-4, 1.0);

  return I;
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(real_t t, N_Vector U, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context)
{
  START_FUNCTION_TIMER();
  am_ode_t* solver = context;
  if (!jacobian_is_current)
  {
    // Compute the approximate Jacobian using a solution vector with ghosts.
    memcpy(solver->U_with_ghosts, NV_DATA(U), sizeof(real_t) * solver->num_local_values);
    newton_pc_setup(solver->precond, 1.0, -gamma, 0.0, t, solver->U_with_ghosts, NULL);
    *jacobian_was_updated = 1;
  }
  else
    *jacobian_was_updated = 0;
  STOP_FUNCTION_TIMER();
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
  START_FUNCTION_TIMER();
  ASSERT(lr == 1); // Left preconditioning only.

  am_ode_t* solver = context;
  
  // FIXME: Apply scaling if needed.

  // Solve it.
  int status = (newton_pc_solve(solver->precond, t, NV_DATA(U), NULL,
                                NV_DATA(r), NV_DATA(r))) ? 0 : 1;
  STOP_FUNCTION_TIMER();
  return status;
}

// Adaptor for J*y setup function
static int set_up_Jy(real_t t, N_Vector y, N_Vector Jy, void* context)
{
  // Nothing here!
  return 0;
}

// Adaptor for J*y evaluation function
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector U, N_Vector U_dot, void* context, N_Vector tmp)
{
  START_FUNCTION_TIMER();
  am_ode_t* solver = context;
  real_t* xx = NV_DATA(U);
  real_t* yy = NV_DATA(y);
  real_t* my_rhs = NV_DATA(U_dot);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);

  // Make sure we use ghosts.
  memcpy(solver->U_with_ghosts, xx, sizeof(real_t) * solver->num_local_values);
  int status = solver->Jy(solver->context, t, solver->U_with_ghosts, my_rhs, yy, temp, jjy);

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < solver->observers->size; ++i)
  {
    am_ode_observer_t* obs = solver->observers->data[i];
    if (obs->Jy_computed != NULL)
      obs->Jy_computed(obs->context, t, xx, my_rhs, yy, jjy);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

ode_solver_t* jfnk_am_ode_solver_new(int order,
                                     MPI_Comm comm,
                                     int num_local_values, 
                                     int num_remote_values, 
                                     void* context, 
                                     int (*rhs_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                                     int (*Jy_func)(void* context, real_t t, real_t* x, real_t* U_dot, real_t* y, real_t* temp, real_t* Jy),
                                     void (*dtor)(void* context),
                                     newton_pc_t* precond,
                                     jfnk_am_krylov_t solver_type,
                                     int max_krylov_dim)
{
  ASSERT(order >= 1);
  ASSERT(order <= 12);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs_func != NULL);
  ASSERT(precond != NULL);
  ASSERT(max_krylov_dim > 3);
  ASSERT(!newton_pc_coefficients_fixed(precond));

  am_ode_t* solver = polymec_malloc(sizeof(am_ode_t));
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

  // Set up KINSol and accessories.
  solver->U = N_VNew(solver->comm, solver->num_local_values);
  solver->U_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->cvode = CVodeCreate(CV_ADAMS);
  CVodeSetMaxOrd(solver->cvode, order);
  CVodeSetUserData(solver->cvode, solver);
  CVodeInit(solver->cvode, am_evaluate_rhs, 0.0, solver->U);
  SUNNonlinearSolver nls = SUNNonlinSol_Newton(solver->U);
  CVodeSetNonlinearSolver(solver->cvode, nls);

  newton_pc_side_t newton_side = newton_pc_side(solver->precond);
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
  if (solver_type == JFNK_AM_GMRES)
  {
    SUNLinearSolver ls = SUNSPGMR(solver->U, side, solver->max_krylov_dim);
    // We use modified Gram-Schmidt orthogonalization.
    SUNSPGMRSetGSType(ls, MODIFIED_GS);
    CVodeSetLinearSolver(solver->cvode, ls, NULL);
  }
  else if (solver_type == JFNK_AM_BICGSTAB)
  {
    SUNLinearSolver ls = SUNSPBCGS(solver->U, side, solver->max_krylov_dim);
    CVodeSetLinearSolver(solver->cvode, ls, NULL);
  }
  else
  {
    SUNLinearSolver ls = SUNSPTFQMR(solver->U, side, solver->max_krylov_dim);
    CVodeSetLinearSolver(solver->cvode, ls, NULL);
  }

  // Set up the Jacobian function and preconditioner.
  if (Jy_func != NULL)
    CVodeSetJacTimes(solver->cvode, set_up_Jy, eval_Jy);
  solver->precond = precond;
  CVodeSetPreconditioner(solver->cvode, set_up_preconditioner,
                           solve_preconditioner_system);

  ode_solver_vtable vtable = {.step = am_step, .advance = am_advance, .dtor = am_dtor};
  char name[1024];
  snprintf(name, 1024, "JFNK Adams-Moulton (order %d)", order);
  ode_solver_t* I = ode_solver_new(name, solver, vtable, order,
                                   num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  am_ode_solver_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void* am_ode_solver_context(ode_solver_t* solver)
{
  am_ode_t* am = ode_solver_context(solver);
  return am->context;
}

void am_ode_solver_set_max_err_test_failures(ode_solver_t* solver,
                                             int max_failures)
{
  ASSERT(max_failures > 0);
  am_ode_t* am = ode_solver_context(solver);
  CVodeSetMaxErrTestFails(am->cvode, max_failures);
}

void am_ode_solver_set_max_nonlinear_iterations(ode_solver_t* solver,
                                                int max_iterations)
{
  ASSERT(max_iterations > 0);
  am_ode_t* am = ode_solver_context(solver);
  CVodeSetMaxNonlinIters(am->cvode, max_iterations);
}

void am_ode_solver_set_nonlinear_convergence_coeff(ode_solver_t* solver,
                                                   real_t coefficient)
{
  ASSERT(coefficient > 0.0);
  am_ode_t* am = ode_solver_context(solver);
  CVodeSetNonlinConvCoef(am->cvode, (double)coefficient);
}

void am_ode_solver_set_tolerances(ode_solver_t* solver,
                                  real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  am_ode_t* am = ode_solver_context(solver);

  // Clear any existing error weight function.
  am->compute_weights = NULL;
  if (am->error_weights != NULL)
  {
    polymec_free(am->error_weights);
    am->error_weights = NULL;
  }

  // Set the tolerances.
  CVodeSStolerances(am->cvode, relative_tol, absolute_tol);
}

// Constant error weight adaptor function.
static void use_constant_weights(void* context, real_t* y, real_t* weights)
{
  am_ode_t* solver = context;
  ASSERT(solver->error_weights != NULL);
  memcpy(weights, solver->error_weights, sizeof(real_t) * solver->num_local_values);
}

void am_ode_solver_set_error_weights(ode_solver_t* solver, real_t* weights)
{
  am_ode_t* am = ode_solver_context(solver);
#ifndef NDEBUG
  // Check for non-negativity and total positivity.
  real_t total = 0.0;
  for (int i = 0; i < am->num_local_values; ++i)
  {
    ASSERT(weights[i] >= 0.0);
    total += weights[i];
  }
  ASSERT(total > 0.0);
#endif

  if (am->error_weights == NULL)
    am->error_weights = polymec_malloc(sizeof(real_t) * am->num_local_values);
  memcpy(am->error_weights, weights, sizeof(real_t) * am->num_local_values);
  am_ode_solver_set_error_weight_function(solver, use_constant_weights);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  am_ode_t* am = context;
  am->compute_weights(am->context, NV_DATA(y), NV_DATA(ewt));

  // Check that all the weights are non-negative.
  int N = (int)NV_LOCLENGTH(y);
  for (int i = 0; i < N; ++i)
  {
    if (NV_Ith(y, i) < 0.0)
      return -1;
  }
  return 0;
}

void am_ode_solver_set_error_weight_function(ode_solver_t* solver,
                                             void (*compute_weights)(void* context, real_t* y, real_t* weights))
{
  am_ode_t* am = ode_solver_context(solver);
  ASSERT(compute_weights != NULL);
  am->compute_weights = compute_weights;
  CVodeWFtolerances(am->cvode, compute_error_weights);
}

void am_ode_solver_eval_rhs(ode_solver_t* solver, real_t t, real_t* U, real_t* U_dot)
{
  am_ode_t* am = ode_solver_context(solver);
  am->rhs(am->context, t, U, U_dot);
}

newton_pc_t* am_ode_solver_preconditioner(ode_solver_t* solver)
{
  am_ode_t* am = ode_solver_context(solver);
  return am->precond;
}

void am_ode_solver_get_diagnostics(ode_solver_t* solver, 
                                   am_ode_solver_diagnostics_t* diagnostics)
{
  am_ode_t* am = ode_solver_context(solver);
  diagnostics->status_message = am->status_message; // borrowed!
  CVodeGetNumSteps(am->cvode, &diagnostics->num_steps);
  CVodeGetLastOrder(am->cvode, &diagnostics->order_of_last_step);
  CVodeGetCurrentOrder(am->cvode, &diagnostics->order_of_next_step);
  CVodeGetLastStep(am->cvode, &diagnostics->last_step_size);
  CVodeGetCurrentStep(am->cvode, &diagnostics->next_step_size);
  CVodeGetNumRhsEvals(am->cvode, &diagnostics->num_rhs_evaluations);
  CVodeGetNumErrTestFails(am->cvode, &diagnostics->num_error_test_failures);
  CVodeGetNumNonlinSolvIters(am->cvode, &diagnostics->num_nonlinear_solve_iterations);
  CVodeGetNumNonlinSolvConvFails(am->cvode, &diagnostics->num_nonlinear_solve_convergence_failures);
  if (am->max_krylov_dim > 0)
  {
    CVodeGetNumLinSolvSetups(am->cvode, &diagnostics->num_linear_solve_setups);
    CVodeGetNumLinIters(am->cvode, &diagnostics->num_linear_solve_iterations);
    CVodeGetNumPrecEvals(am->cvode, &diagnostics->num_preconditioner_evaluations);
    CVodeGetNumPrecSolves(am->cvode, &diagnostics->num_preconditioner_solves);
    CVodeGetNumLinConvFails(am->cvode, &diagnostics->num_linear_solve_convergence_failures);
  }
  else 
  {
    diagnostics->num_linear_solve_setups = -1;
    diagnostics->num_linear_solve_iterations = -1;
    diagnostics->num_preconditioner_evaluations = -1;
    diagnostics->num_preconditioner_solves = -1;
    diagnostics->num_linear_solve_convergence_failures = -1;
  }
}

void am_ode_solver_diagnostics_fprintf(am_ode_solver_diagnostics_t* diagnostics, 
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
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  if (diagnostics->num_linear_solve_setups != -1)
    fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  if (diagnostics->num_linear_solve_iterations != -1)
    fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  if (diagnostics->num_linear_solve_convergence_failures != -1)
    fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  if (diagnostics->num_preconditioner_evaluations != -1)
    fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  if (diagnostics->num_preconditioner_solves != -1)
    fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
}

am_ode_observer_t* am_ode_observer_new(void* context,
                                       void (*rhs_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                       void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy),
                                       void (*dtor)(void* context))
{
  am_ode_observer_t* obs = polymec_malloc(sizeof(am_ode_observer_t));
  obs->context = context;
  obs->rhs_computed = rhs_computed;
  obs->Jy_computed = Jy_computed;
  obs->dtor = dtor;
  return obs;
}

static void am_ode_observer_free(am_ode_observer_t* observer)
{
  if ((observer->dtor != NULL) && (observer->context != NULL))
    observer->dtor(observer->context);
  polymec_free(observer);
}

void am_ode_solver_add_observer(ode_solver_t* solver,
                                am_ode_observer_t* observer)
{
  am_ode_t* am = ode_solver_context(solver);
  ptr_array_append_with_dtor(am->observers, observer, DTOR(am_ode_observer_free));
}

