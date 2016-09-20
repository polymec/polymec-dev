// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "core/array.h"
#include "core/timer.h"
#include "integrators/bdf_ode_integrator.h"
#include "integrators/newton_pc.h"
#include "cvode/cvode.h"

// JFNK stuff.
#include "cvode/cvode_spils.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"
#include "cvode/cvode_sptfqmr.h"

// Stuff for generalized BDF integrators.
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

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  int (*Jy)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* temp, real_t* Jy);

  // Generalized adaptor stuff.
  int (*reset_func)(void* context, real_t t, real_t* X);
  int (*setup_func)(void* context, 
                    bdf_conv_status_t conv_status, 
                    real_t gamma, 
                    real_t t, 
                    real_t* U_pred, 
                    real_t* U_dot_pred, 
                    bool* J_current, 
                    real_t* work1, real_t* work2, real_t* work3);
  int (*solve_func)(void* context, 
                    real_t* W, 
                    real_t t, 
                    real_t* U,
                    real_t* U_dot,
                    real_t* B);

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
    snprintf(err, 1024, "Integrator stopped at t = %g after maximum number of steps.", current_time);
    status_message = string_dup(err);
  }
  else if (status == CV_TOO_MUCH_ACC)
    status_message = string_dup("Integrator could not achieve desired level of accuracy.");
  else if (status == CV_ERR_FAILURE)
    status_message = string_dup("Integrator encountered too many error test failures.");
  else if (status == CV_CONV_FAILURE)
    status_message = string_dup("Integrator encountered too many convergence test failures.");
  else if (status == CV_LINIT_FAIL)
    status_message = string_dup("Integrator's linear solver failed to initialize.");
  else if (status == CV_LSETUP_FAIL)
    status_message = string_dup("Integrator's linear solver setup failed.");
  else if (status == CV_LSOLVE_FAIL)
    status_message = string_dup("Integrator's linear solver failed.");
  else if (status == CV_RHSFUNC_FAIL)
    status_message = string_dup("Integrator's RHS function failed unrecoverably.");
  else if (status == CV_FIRST_RHSFUNC_ERR)
    status_message = string_dup("Integrator's first call to RHS function failed.");
  else if (status == CV_REPTD_RHSFUNC_ERR)
    status_message = string_dup("Integrator encountered too many recoverable RHS failures.");
  else if (status == CV_UNREC_RHSFUNC_ERR)
    status_message = string_dup("Integrator failed to recover from a recoverable RHS failure.");
  else if (status == CV_RTFUNC_FAIL)
    status_message = string_dup("Integrator encountered a failure in the rootfinding function.");
  return status_message;
}

// This function wraps around the user-supplied right hand side.
static int bdf_evaluate_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  START_FUNCTION_TIMER();
  bdf_ode_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* xxd = NV_DATA(U_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->rhs(integ->context, t, integ->U_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, integ->U_with_ghosts, sizeof(real_t) * integ->num_local_values);
  }

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < integ->observers->size; ++i)
  {
    bdf_ode_observer_t* obs = integ->observers->data[i];
    if (obs->rhs_computed != NULL)
      obs->rhs_computed(obs->context, t, xx, xxd);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

static bool bdf_step(void* context, real_t max_dt, real_t* t, real_t* U)
{
  START_FUNCTION_TIMER();
  bdf_ode_t* integ = context;
  int status = CV_SUCCESS;

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > integ->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = CVode(integ->cvode, t2, integ->U, &integ->t, CV_ONE_STEP);
    if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN))
    {
      integ->status_message = get_status_message(status, integ->t);
      STOP_FUNCTION_TIMER();
      return false;
    }

    if ((t2 - *t) < (integ->t - *t))
      log_detail("bdf_ode_integrator: took internal step dt = %g", integ->t - *t);
  }

  // If we integrated past t2, interpolate to t2.
  if (integ->t > t2)
  {
    status = CVodeGetDky(integ->cvode, t2, 0, integ->U);
    *t = t2;
  }
  else
    *t = integ->t;
  
  // Clear the present status.
  if (integ->status_message != NULL)
  {
    polymec_free(integ->status_message);
    integ->status_message = NULL;
  }

  // Did it work?
  if ((status == CV_SUCCESS) || (status == CV_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    STOP_FUNCTION_TIMER();
    return false;
  }
}

static bool bdf_advance(void* context, real_t t1, real_t t2, real_t* U)
{
  bdf_ode_t* integ = context;
  integ->t = t1;
  CVodeReInit(integ->cvode, t1, integ->U);
  CVodeSetStopTime(integ->cvode, t2);
  
  // Copy in the solution.
  memcpy(NV_DATA(integ->U), U, sizeof(real_t) * integ->num_local_values); 

  // Integrate.
  int status = CVode(integ->cvode, t2, integ->U, &integ->t, CV_NORMAL);
  
  // Clear the present status.
  if (integ->status_message != NULL)
  {
    polymec_free(integ->status_message);
    integ->status_message = NULL;
  }

  // Did it work?
  if ((status == CV_SUCCESS) || (status == CV_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static void bdf_reset(void* context, real_t t, real_t* U)
{
  bdf_ode_t* integ = context;

  // Reset the preconditioner.
  if (integ->precond != NULL)
    newton_pc_reset(integ->precond, t);

  // Copy in the solution and reinitialize.
  memcpy(NV_DATA(integ->U), U, sizeof(real_t) * integ->num_local_values); 
  CVodeReInit(integ->cvode, t, integ->U);
  integ->t = t;
}

static void bdf_dtor(void* context)
{
  bdf_ode_t* integ = context;

  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the CVode stuff.
  polymec_free(integ->U_with_ghosts);
  N_VDestroy(integ->U);
  CVodeFree(&integ->cvode);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  ptr_array_free(integ->observers);
  if (integ->error_weights != NULL)
    polymec_free(integ->error_weights);
  polymec_free(integ);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector U, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  bdf_ode_t* integ = context;
  if (!jacobian_is_current)
  {
    // Compute the approximate Jacobian using a solution vector with ghosts.
    memcpy(integ->U_with_ghosts, NV_DATA(U), sizeof(real_t) * integ->num_local_values);
    newton_pc_setup(integ->precond, 1.0, -gamma, 0.0, t, integ->U_with_ghosts, NULL);
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
                                       int lr, void* context, 
                                       N_Vector work)
{
  bdf_ode_t* integ = context;
  
  // FIXME: Apply scaling if needed.

  // Set the preconditioner's L2 norm tolerance.
  newton_pc_set_tolerance(integ->precond, delta);

  // Solve it.
  int result = (newton_pc_solve(integ->precond, t, NV_DATA(U), NULL,
                                NV_DATA(r), NV_DATA(z))) ? 0 : 1;
  return result;
}

// Adaptor for J*y function
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector U, N_Vector rhs, void* context, N_Vector tmp)
{
  START_FUNCTION_TIMER();
  bdf_ode_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* yy = NV_DATA(y);
  real_t* my_rhs = NV_DATA(rhs);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);

  // Make sure we use ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->Jy(integ->context, t, xx, my_rhs, yy, temp, jjy);

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < integ->observers->size; ++i)
  {
    bdf_ode_observer_t* obs = integ->observers->data[i];
    if (obs->Jy_computed != NULL)
      obs->Jy_computed(obs->context, t, xx, my_rhs, yy, jjy);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

ode_integrator_t* jfnk_bdf_ode_integrator_new(int order,
                                              MPI_Comm comm,
                                              int num_local_values, 
                                              int num_remote_values, 
                                              void* context, 
                                              int (*rhs_func)(void* context, real_t t, real_t* U, real_t* xdot),
                                              int (*Jy_func)(void* context, real_t t, real_t* U, real_t* rhs, real_t* y, real_t* temp, real_t* Jy),
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

  bdf_ode_t* integ = polymec_malloc(sizeof(bdf_ode_t));
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->rhs = rhs_func;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = max_krylov_dim;
  integ->Jy = Jy_func;
  integ->t = 0.0;
  integ->observers = ptr_array_new();
  integ->error_weights = NULL;

  integ->reset_func = NULL;
  integ->setup_func = NULL;
  integ->solve_func = NULL;

  // Set up CVode and accessories.
  integ->U = N_VNew(integ->comm, integ->num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
  integ->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetMaxOrd(integ->cvode, order);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, bdf_evaluate_rhs, 0.0, integ->U);

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
    CVSpgmr(integ->cvode, side, max_krylov_dim); 
    // We use modified Gram-Schmidt orthogonalization.
    CVSpilsSetGSType(integ->cvode, MODIFIED_GS);
  }
  else if (solver_type == JFNK_BDF_BICGSTAB)
    CVSpbcg(integ->cvode, side, max_krylov_dim);
  else
    CVSptfqmr(integ->cvode, side, max_krylov_dim);

  // Set up the Jacobian function and preconditioner.
  if (Jy_func != NULL)
    CVSpilsSetJacTimesVecFn(integ->cvode, eval_Jy);
  integ->precond = precond;
  CVSpilsSetPreconditioner(integ->cvode, set_up_preconditioner,
                           solve_preconditioner_system);

  ode_integrator_vtable vtable = {.step = bdf_step, .advance = bdf_advance, .reset = bdf_reset, .dtor = bdf_dtor};
  char name[1024];
  snprintf(name, 1024, "JFNK Backwards-Difference-Formulae (order %d)", order);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order,
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void* bdf_ode_integrator_context(ode_integrator_t* integrator)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  return integ->context;
}

void bdf_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                  int max_failures)
{
  ASSERT(max_failures > 0);
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetMaxErrTestFails(integ->cvode, max_failures);
}

void bdf_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                     int max_iterations)
{
  ASSERT(max_iterations > 0);
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetMaxNonlinIters(integ->cvode, max_iterations);
}

void bdf_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                        real_t coefficient)
{
  ASSERT(coefficient > 0.0);
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetNonlinConvCoef(integ->cvode, (double)coefficient);
}

void bdf_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                      real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  bdf_ode_t* integ = ode_integrator_context(integrator);

  // Clear any existing error weight function.
  integ->compute_weights = NULL;
  if (integ->error_weights != NULL)
  {
    polymec_free(integ->error_weights);
    integ->error_weights = NULL;
  }

  // Set the tolerances.
  CVodeSStolerances(integ->cvode, relative_tol, absolute_tol);
}

// Constant error weight adaptor function.
static void use_constant_weights(void* context, real_t* y, real_t* weights)
{
  bdf_ode_t* integ = context;
  ASSERT(integ->error_weights != NULL);
  memcpy(weights, integ->error_weights, sizeof(real_t) * integ->num_local_values);
}

void bdf_ode_integrator_set_error_weights(ode_integrator_t* integrator, real_t* weights)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
#ifndef NDEBUG
  // Check for non-negativity and total positivity.
  real_t total = 0.0;
  for (int i = 0; i < integ->num_local_values; ++i)
  {
    ASSERT(weights[i] >= 0.0);
    total += weights[i];
  }
  ASSERT(total > 0.0);
#endif

  if (integ->error_weights == NULL)
    integ->error_weights = polymec_malloc(sizeof(real_t) * integ->num_local_values);
  memcpy(integ->error_weights, weights, sizeof(real_t) * integ->num_local_values);
  bdf_ode_integrator_set_error_weight_function(integrator, use_constant_weights);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  bdf_ode_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));

  // Check that all the weights are non-negative.
  int N = (int)NV_LOCLENGTH(y);
  for (int i = 0; i < N; ++i)
  {
    if (NV_Ith(y, i) < 0.0)
      return -1;
  }
  return 0;
}

void bdf_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                  void (*compute_weights)(void* context, real_t* y, real_t* weights))
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  ASSERT(compute_weights != NULL);
  integ->compute_weights = compute_weights;
  CVodeWFtolerances(integ->cvode, compute_error_weights);
}

void bdf_ode_integrator_set_stability_limit_detection(ode_integrator_t* integrator,
                                                      bool use_detection)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetStabLimDet(integ->cvode, use_detection);
}

void bdf_ode_integrator_eval_rhs(ode_integrator_t* integrator, real_t t, real_t* U, real_t* U_dot)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  memcpy(integ->U_with_ghosts, U, sizeof(real_t) * integ->num_local_values);
  integ->rhs(integ->context, t, integ->U_with_ghosts, U_dot);
}

newton_pc_t* bdf_ode_integrator_preconditioner(ode_integrator_t* integrator)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  return integ->precond;
}

void bdf_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                        bdf_ode_integrator_diagnostics_t* diagnostics)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  diagnostics->status_message = integ->status_message; // borrowed!
  CVodeGetNumSteps(integ->cvode, &diagnostics->num_steps);
  CVodeGetLastOrder(integ->cvode, &diagnostics->order_of_last_step);
  CVodeGetCurrentOrder(integ->cvode, &diagnostics->order_of_next_step);
  CVodeGetLastStep(integ->cvode, &diagnostics->last_step_size);
  CVodeGetCurrentStep(integ->cvode, &diagnostics->next_step_size);
  CVodeGetNumRhsEvals(integ->cvode, &diagnostics->num_rhs_evaluations);
  CVodeGetNumLinSolvSetups(integ->cvode, &diagnostics->num_linear_solve_setups);
  CVodeGetNumErrTestFails(integ->cvode, &diagnostics->num_error_test_failures);
  CVodeGetNumNonlinSolvIters(integ->cvode, &diagnostics->num_nonlinear_solve_iterations);
  CVodeGetNumNonlinSolvConvFails(integ->cvode, &diagnostics->num_nonlinear_solve_convergence_failures);
  if (integ->solve_func == NULL) // JFNK mode
  {
    CVSpilsGetNumLinIters(integ->cvode, &diagnostics->num_linear_solve_iterations);
    CVSpilsGetNumPrecEvals(integ->cvode, &diagnostics->num_preconditioner_evaluations);
    CVSpilsGetNumPrecSolves(integ->cvode, &diagnostics->num_preconditioner_solves);
    CVSpilsGetNumConvFails(integ->cvode, &diagnostics->num_linear_solve_convergence_failures);
  }
  else
  {
    diagnostics->num_linear_solve_iterations = -1;
    diagnostics->num_preconditioner_evaluations = -1;
    diagnostics->num_preconditioner_solves = -1;
    diagnostics->num_linear_solve_convergence_failures = -1;
  }
}

void bdf_ode_integrator_diagnostics_fprintf(bdf_ode_integrator_diagnostics_t* diagnostics, 
                                            FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "ODE integrator diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Order of last step: %d\n", diagnostics->order_of_last_step);
  fprintf(stream, "  Order of next step: %d\n", diagnostics->order_of_next_step);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Next step size: %g\n", diagnostics->next_step_size);
  fprintf(stream, "  Num RHS evaluations: %d\n", (int)diagnostics->num_rhs_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  if (diagnostics->num_linear_solve_convergence_failures != -1) // JFNK mode
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

void bdf_ode_integrator_add_observer(ode_integrator_t* integrator,
                                     bdf_ode_observer_t* observer)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  ptr_array_append_with_dtor(integ->observers, observer, DTOR(bdf_ode_observer_free));
}

// This unpublished function provides direct access to the CVode object 
// within the BDF integrator.
void* bdf_ode_integrator_cvode(ode_integrator_t* integrator);
void* bdf_ode_integrator_cvode(ode_integrator_t* integrator)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  return integ->cvode;
}

// This unpublished function provides direct access to the CVode right-hand
// side function within the BDF integrator.
int bdf_ode_integrator_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context);
int bdf_ode_integrator_rhs(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  return bdf_evaluate_rhs(t, U, U_dot, context);
}

//------------------------------------------------------------------------
//                    Custom BDF integrator stuff
//------------------------------------------------------------------------

static int bdf_linit(CVodeMem cv_mem)
{
  bdf_ode_t* bdf = cv_mem->cv_user_data;
  real_t t = cv_mem->cv_tn;
  real_t* U = NV_DATA(cv_mem->cv_y);
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
  {
    if (cv_mem->cv_nst == 0)
      conv_status = BDF_CONV_FIRST_STEP;
    else
      conv_status = BDF_CONV_NO_FAILURES;
  }
  else if (convfail == CV_FAIL_BAD_J)
    conv_status = BDF_CONV_BAD_J_FAILURE;
  else
    conv_status = BDF_CONV_OTHER_FAILURE;
  real_t gamma = cv_mem->cv_gamma;
  real_t t = cv_mem->cv_tn;
  bool J_current = false;
  real_t* U_pred = NV_DATA(ypred);
  real_t* U_dot_pred = NV_DATA(fpred);
  real_t* work1 = NV_DATA(vtemp1);
  real_t* work2 = NV_DATA(vtemp2);
  real_t* work3 = NV_DATA(vtemp3);
  int status = bdf->setup_func(bdf->context, conv_status, gamma, t, 
                               U_pred, U_dot_pred, &J_current, work1, 
                               work2, work3);
  *jcurPtr = J_current;
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
  real_t* W = NV_DATA(weight);
  real_t* U = NV_DATA(ycur);
  real_t* U_dot = NV_DATA(fcur);
  real_t* B = NV_DATA(b);
  return bdf->solve_func(bdf->context, W, t, U, U_dot, B);
}

static void bdf_lfree(CVodeMem cv_mem)
{
}

ode_integrator_t* bdf_ode_integrator_new(const char* name,
                                         int order, 
                                         MPI_Comm comm,
                                         int num_local_values, 
                                         int num_remote_values, 
                                         void* context, 
                                         int (*rhs_func)(void* context, real_t t, real_t* U, real_t* xdot),
                                         int (*reset_func)(void* context, real_t t, real_t* X),
                                         int (*setup_func)(void* context, 
                                                           bdf_conv_status_t conv_status, 
                                                           real_t gamma, 
                                                           real_t t, 
                                                           real_t* X_pred, 
                                                           real_t* rhs_pred, 
                                                           bool* J_current, 
                                                           real_t* work1, real_t* work2, real_t* work3),
                                         int (*solve_func)(void* context, 
                                                           real_t* W, 
                                                           real_t t, 
                                                           real_t* X,
                                                           real_t* rhs,
                                                           real_t* B),
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

  bdf_ode_t* integ = polymec_malloc(sizeof(bdf_ode_t));
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->rhs = rhs_func;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = -1;
  integ->Jy = NULL;
  integ->precond = NULL;
  integ->t = 0.0;
  integ->observers = ptr_array_new();
  integ->error_weights = NULL;

  // Set up CVode and accessories.
  integ->U = N_VNew(integ->comm, integ->num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
  integ->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetMaxOrd(integ->cvode, order);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, bdf_evaluate_rhs, 0.0, integ->U);

  // Set up the solver.
  integ->reset_func = reset_func;
  integ->setup_func = setup_func;
  integ->solve_func = solve_func;
  {
    CVodeMem cv_mem = integ->cvode;
    cv_mem->cv_linit = bdf_linit;
    cv_mem->cv_lsetup = bdf_lsetup;
    cv_mem->cv_lsolve = bdf_lsolve;
    cv_mem->cv_lfree = bdf_lfree;
    cv_mem->cv_setupNonNull = true;
  }

  ode_integrator_vtable vtable = {.step = bdf_step, 
                                  .advance = bdf_advance, 
                                  .reset = bdf_reset, 
                                  .dtor = bdf_dtor};
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order,
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

//------------------------------------------------------------------------
//                  Inexact Newton-Krylov integrator stuff
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;

  void* context;
  int (*J_func)(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J);
  void (*dtor)(void* context);

  krylov_factory_t* factory;
  krylov_solver_t* solver;
  krylov_pc_t* pc;
  krylov_matrix_t* J;
  krylov_vector_t* X;
  krylov_vector_t* B;

  matrix_sparsity_t* sparsity;
  int block_size;
  index_t start, end;
} ink_bdf_ode_t;

static int ink_reset(void* context, real_t t, real_t* X)
{
  START_FUNCTION_TIMER();
  // Allocate resources.
  ink_bdf_ode_t* ink = context;
  if (ink->J != NULL)
    krylov_matrix_free(ink->J);
  if (ink->block_size > 1)
    ink->J = krylov_factory_block_matrix(ink->factory, ink->sparsity, ink->block_size);
  else
    ink->J = krylov_factory_matrix(ink->factory, ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);
  index_t* row_dist = matrix_sparsity_row_distribution(ink->sparsity);
  ink->X = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->B = krylov_factory_vector(ink->factory, ink->comm, row_dist);

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_setup(void* context, 
                     bdf_conv_status_t conv_status, 
                     real_t gamma, 
                     real_t t, 
                     real_t* U_pred, 
                     real_t* U_dot_pred, 
                     bool* J_updated, 
                     real_t* work1, real_t* work2, real_t* work3)
{
  START_FUNCTION_TIMER();
  ink_bdf_ode_t* ink = context;

  if ((conv_status == BDF_CONV_FIRST_STEP) || 
      (conv_status == BDF_CONV_BAD_J_FAILURE) || 
      (conv_status == BDF_CONV_OTHER_FAILURE))
  {
    // Call our Jacobian calculation function.
    log_debug("ink_bdf_ode_integrator: Calculating A = I - %g * J.\n", gamma);
    int status = ink->J_func(ink->context, t, U_pred, U_dot_pred, ink->J);
    if (status != 0)
      return status;

    // Scale the Jacobian by -gamma and add the identity matrix.
    krylov_matrix_scale(ink->J, -gamma);
    krylov_matrix_add_identity(ink->J, 1.0);

    // Use this matrix as the operator in our solver.
    krylov_solver_set_operator(ink->solver, ink->J);

    *J_updated = true;
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    // The previous Newton iteration failed to converge even with a current Jacobian.
    *J_updated = false;
    STOP_FUNCTION_TIMER();
    return (conv_status == BDF_CONV_NO_FAILURES) ? 0 : -1;
  }
}

static int ink_solve(void* context, 
                     real_t* W, 
                     real_t t, 
                     real_t* U,
                     real_t* U_dot,
                     real_t* B) 
{
  START_FUNCTION_TIMER();
  ink_bdf_ode_t* ink = context;

  // Copy RHS data from B into ink->B.
  krylov_vector_copy_in(ink->B, B);

  real_t res_norm;
  int num_iters;
  bool solved = krylov_solver_solve(ink->solver, ink->B, ink->X, &res_norm, &num_iters);

  if (solved)
  {
    log_debug("ink_bdf_ode_integrator: Solved A*X = B (||R|| == %g after %d iters).", res_norm, num_iters);

    // Copy solution data from ink->X into B.
    krylov_vector_copy_out(ink->X, B);
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    log_debug("ink_bdf_ode_integrator: Solution to A*X = B did not converge.");
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
  if (ink->J != NULL)
    krylov_matrix_free(ink->J);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  krylov_solver_free(ink->solver);
  krylov_factory_free(ink->factory);
  if ((ink->context != NULL) && (ink->dtor != NULL))
    ink->dtor(ink->context);
}

ode_integrator_t* ink_bdf_ode_integrator_new(int order, 
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
  ink->J_func = J_func;
  ink->dtor = dtor;
  ink->factory = factory;
  ink->sparsity = J_sparsity;
  ink->solver = krylov_factory_gmres_solver(ink->factory, comm, 30);
  ink->pc = NULL;
  ink->J = NULL;
  ink->B = NULL;
  ink->X = NULL;

  // Find the start, end indices.
  int rank;
  MPI_Comm_rank(comm, &rank);
  index_t* row_dist = matrix_sparsity_row_distribution(J_sparsity);
  ink->start = row_dist[rank];
  ink->end = row_dist[rank+1];
  ink->block_size = 1;

  char name[1024];
  snprintf(name, 1024, "INK Backwards-Difference-Formulae (order %d)", order);
  int num_local_values = (int)(matrix_sparsity_num_local_rows(J_sparsity));
  ode_integrator_t* I = bdf_ode_integrator_new(name, order, comm, 
                                               num_local_values, 0,
                                               ink, rhs_func, ink_reset, 
                                               ink_setup, ink_solve, ink_dtor);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void ink_bdf_ode_integrator_use_pcg(ode_integrator_t* ink_bdf_ode_integ)
{
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_pcg_solver(ink->factory, ink->comm);
}

void ink_bdf_ode_integrator_use_gmres(ode_integrator_t* ink_bdf_ode_integ,
                                      int max_krylov_dim)
{
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_gmres_solver(ink->factory, ink->comm, max_krylov_dim);
}

void ink_bdf_ode_integrator_use_bicgstab(ode_integrator_t* ink_bdf_ode_integ)
{
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_bicgstab_solver(ink->factory, ink->comm);
}

void ink_bdf_ode_integrator_use_special(ode_integrator_t* ink_bdf_ode_integ,
                                        const char* solver_name,
                                        string_string_unordered_map_t* options)
{
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_special_solver(ink->factory, ink->comm,
                                              solver_name, options);
}

void ink_bdf_ode_integrator_set_pc(ode_integrator_t* ink_bdf_ode_integ,
                                   const char* pc_name, 
                                   string_string_unordered_map_t* options)
{
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  ink->pc = krylov_factory_preconditioner(ink->factory, ink->comm, pc_name, options);
}

void ink_bdf_ode_integrator_set_block_size(ode_integrator_t* ink_bdf_ode_integ,
                                           int block_size)
{
  ASSERT(block_size > 0);
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  ink->block_size = block_size;
}

void* ink_bdf_ode_integrator_context(ode_integrator_t* ink_bdf_ode_integ)
{
  ink_bdf_ode_t* ink = bdf_ode_integrator_context(ink_bdf_ode_integ);
  return ink->context;
}
