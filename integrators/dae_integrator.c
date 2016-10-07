// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "core/timer.h"
#include "integrators/dae_integrator.h"

#include "ida/ida.h"
#include "ida/ida_spils.h"
#include "ida/ida_spgmr.h"
#include "ida/ida_spbcgs.h"
#include "ida/ida_sptfqmr.h"

// These symbols represent special short-hand arguments for equation types 
// and constraints.
static dae_equation_t DAE_ALL_ALGEBRAIC_LOC = DAE_ALGEBRAIC;
dae_equation_t* DAE_ALL_ALGEBRAIC = &DAE_ALL_ALGEBRAIC_LOC;

static dae_equation_t DAE_ALL_DIFFERENTIAL_LOC = DAE_DIFFERENTIAL;
dae_equation_t* DAE_ALL_DIFFERENTIAL = &DAE_ALL_DIFFERENTIAL_LOC;

static dae_constraint_t DAE_ALL_UNCONSTRAINED_LOC = DAE_UNCONSTRAINED;
dae_constraint_t* DAE_ALL_UNCONSTRAINED = &DAE_ALL_UNCONSTRAINED_LOC;

static dae_constraint_t DAE_ALL_NEGATIVE_LOC = DAE_NEGATIVE;
dae_constraint_t* DAE_ALL_NEGATIVE = &DAE_ALL_NEGATIVE_LOC;

static dae_constraint_t DAE_ALL_NONPOSITIVE_LOC = DAE_NONPOSITIVE;
dae_constraint_t* DAE_ALL_NONPOSITIVE = &DAE_ALL_NONPOSITIVE_LOC;

static dae_constraint_t DAE_ALL_NONNEGATIVE_LOC = DAE_NONNEGATIVE;
dae_constraint_t* DAE_ALL_NONNEGATIVE = &DAE_ALL_NONNEGATIVE_LOC;

static dae_constraint_t DAE_ALL_POSITIVE_LOC = DAE_POSITIVE;
dae_constraint_t* DAE_ALL_POSITIVE = &DAE_ALL_POSITIVE_LOC;

struct dae_integrator_t 
{
  char* name;
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context;
  int order;
  bool initialized;
  int (*F)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F);
  void (*dtor)(void* context);

  // IDA data structures.
  void* ida;
  N_Vector U, U_dot;
  real_t* U_with_ghosts;
  real_t* U_dot_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.
  real_t max_dt, stop_time;

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  jfnk_dae_krylov_t solver_type;
  int (*Jy)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F,
            real_t* y, real_t* Jy, real_t* tmp1, real_t* tmp2);

  // Generalized adaptor stuff.
  int (*reset_func)(void* context, real_t t, real_t* U, real_t* U_dot);
  int (*setup_func)(void* context, 
                    real_t alpha, 
                    int step,
                    real_t t, 
                    real_t* U_pred, 
                    real_t* U_dot_pred, 
                    real_t* F_pred, 
                    real_t* work1, real_t* work2, real_t* work3);
  int (*solve_func)(void* context, 
                    real_t t, 
                    real_t* U,
                    real_t* U_dot,
                    real_t* F,
                    real_t* W, 
                    real_t res_norm_tol,
                    real_t* B);
  int prev_J_step;

  // Error weight function.
  void (*compute_weights)(void* context, real_t* y, real_t* weights);
  real_t* error_weights;
};

static char* get_status_message(int status, real_t current_time)
{
  char* status_message = NULL;
  if (status == IDA_TOO_MUCH_WORK)
  {
    char err[1024];
    snprintf(err, 1024, "Integrator stopped at t = %g after maximum number of steps.", current_time);
    status_message = string_dup(err);
  }
  else if (status == IDA_TOO_MUCH_ACC)
    status_message = string_dup("Integrator could not achieve desired level of accuracy.");
  else if (status == IDA_ERR_FAIL)
    status_message = string_dup("Integrator encountered too many error test failures.");
  else if (status == IDA_CONV_FAIL)
    status_message = string_dup("Integrator encountered too many convergence test failures.");
  else if (status == IDA_LINIT_FAIL)
    status_message = string_dup("Integrator's linear solver failed to initialize.");
  else if (status == IDA_LSETUP_FAIL)
    status_message = string_dup("Integrator's linear solver setup failed.");
  else if (status == IDA_LSOLVE_FAIL)
    status_message = string_dup("Integrator's linear solver failed.");
  else if (status == IDA_CONSTR_FAIL)
    status_message = string_dup("Integrator's inequality constraints were violated unrecoverably.");
  else if (status == IDA_RES_FAIL)
    status_message = string_dup("Integrator's residual function failed unrecoverably.");
  else if (status == IDA_REP_RES_ERR)
    status_message = string_dup("Integrator encountered too many recoverable RHS failures.");
  else if (status == IDA_RTFUNC_FAIL)
    status_message = string_dup("Integrator encountered a failure in the rootfinding function.");
  return status_message;
}

// This function wraps around the user-supplied right hand side.
static int evaluate_residual(real_t t, N_Vector U, N_Vector U_dot, 
                             N_Vector F, void* context)
{
  dae_integrator_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* xxd = NV_DATA(U_dot);
  real_t* Fx = NV_DATA(F);

  // Evaluate the residual using vectors with ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, xxd, sizeof(real_t) * integ->num_local_values);
  return integ->F(integ->context, t, integ->U_with_ghosts, integ->U_dot_with_ghosts, Fx);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector U, N_Vector U_dot, N_Vector F,
                                 real_t cj, void* context,
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  START_FUNCTION_TIMER();
  dae_integrator_t* integ = context;

  // Form the linear combination dFdx + cj * dFd(U_dot).
  memcpy(integ->U_with_ghosts, NV_DATA(U), sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, NV_DATA(U_dot), sizeof(real_t) * integ->num_local_values);
  newton_pc_setup(integ->precond, 0.0, 1.0, cj, t, integ->U_with_ghosts, integ->U_dot_with_ghosts);
  STOP_FUNCTION_TIMER();
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector U, N_Vector U_dot, 
                                       N_Vector F, N_Vector r, N_Vector z, 
                                       real_t cj, real_t delta, 
                                       void* context, N_Vector work)
{
  START_FUNCTION_TIMER();
  dae_integrator_t* integ = context;
  
  // FIXME: Apply scaling if needed.

  // Solve it.
  int result = (newton_pc_solve(integ->precond, t, NV_DATA(U), NV_DATA(U_dot),
                                NV_DATA(r), NV_DATA(z))) ? 0 : 1;
  STOP_FUNCTION_TIMER();
  return result;
}

// Adaptor for J*y function.
static int eval_Jy(real_t tt, N_Vector uu, N_Vector up, N_Vector ff,
                   N_Vector y, N_Vector Jy, real_t cj, 
                   void *context, N_Vector tmp1, N_Vector tmp2)
{
  START_FUNCTION_TIMER();
  dae_integrator_t* integ = context;
  real_t* U = NV_DATA(uu);
  real_t* Udot = NV_DATA(up);
  real_t* F = NV_DATA(ff);
  real_t* yy = NV_DATA(y);
  real_t* Jyy = NV_DATA(Jy);
  real_t* tmp11 = NV_DATA(tmp1);
  real_t* tmp22 = NV_DATA(tmp2);

  // Make sure we use ghosts.
  memcpy(integ->U_with_ghosts, U, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, Udot, sizeof(real_t) * integ->num_local_values);
  int status = integ->Jy(integ->context, tt, U, cj, Udot, F, yy, Jyy, tmp11, tmp22);
  STOP_FUNCTION_TIMER();
  return status;
}

// We use floating point comparisons to assign inequality relations in this 
// function, so we disable related compiler warnings.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
dae_integrator_t* jfnk_dae_integrator_new(int order,
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
                                          int max_krylov_dim)
{
  ASSERT(order > 0);
  ASSERT(order <= 5);
  ASSERT(equation_types != NULL);
  ASSERT(constraints != NULL);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(F_func != NULL);
  ASSERT(precond != NULL);
  ASSERT(newton_pc_side(precond) == NEWTON_PC_LEFT); // left PC only!
  ASSERT(max_krylov_dim >= 3);

  dae_integrator_t* integ = polymec_malloc(sizeof(dae_integrator_t));
  integ->context = context;
  integ->comm = comm;
  integ->order = order;
  integ->solver_type = solver_type;
  integ->t = 0.0;
  integ->F = F_func;
  integ->Jy = Jy_func;
  integ->dtor = dtor;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->precond = precond;
  integ->max_krylov_dim = max_krylov_dim;
  integ->initialized = false;
  integ->max_dt = REAL_MAX;
  integ->status_message = NULL;
  integ->error_weights = NULL;

  integ->reset_func = NULL;
  integ->setup_func = NULL;
  integ->solve_func = NULL;

  // Set up IDA and accessories.
  integ->U = N_VNew(comm, num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->U_dot = N_VNew(comm, num_local_values);
  integ->U_dot_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  memset(NV_DATA(integ->U), 0, sizeof(real_t) * num_local_values);
  memset(NV_DATA(integ->U_dot), 0, sizeof(real_t) * num_local_values);

  integ->ida = IDACreate();
  IDASetMaxOrd(integ->ida, integ->order);
  IDASetUserData(integ->ida, integ);
  IDAInit(integ->ida, evaluate_residual, integ->t, integ->U, integ->U_dot);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver_type == JFNK_DAE_GMRES)
    IDASpgmr(integ->ida, max_krylov_dim); 
  else if (solver_type == JFNK_DAE_BICGSTAB)
    IDASpbcg(integ->ida, max_krylov_dim);
  else
    IDASptfqmr(integ->ida, max_krylov_dim);

  // We set the equation types.
  N_Vector array = N_VNew(comm, num_local_values);
  bool has_algebraic_eqns = false;
  if (equation_types == DAE_ALL_ALGEBRAIC)
  {
    has_algebraic_eqns = true;
    for (int i = 0; i < num_local_values; ++i)
      NV_Ith(array, i) = 0.0;
  }
  else if (equation_types == DAE_ALL_DIFFERENTIAL)
  {
    for (int i = 0; i < num_local_values; ++i)
      NV_Ith(array, i) = 1.0;
  }
  else
  {
    for (int i = 0; i < num_local_values; ++i)
    {
      if (equation_types[i] == DAE_ALGEBRAIC)
      {
        has_algebraic_eqns = true;
        NV_Ith(array, i) = 0.0;
      }
      else
        NV_Ith(array, i) = 1.0;
    }
  }
  IDASetId(integ->ida, array);

  // Set up constraints.
  if (constraints != DAE_ALL_UNCONSTRAINED)
  {
    if (constraints == DAE_ALL_NEGATIVE)
    {
      for (int i = 0; i < num_local_values; ++i)
        NV_Ith(array, i) = -2.0;
    }
    else if (constraints == DAE_ALL_NONPOSITIVE)
    {
      for (int i = 0; i < num_local_values; ++i)
        NV_Ith(array, i) = -1.0;
    }
    else if (constraints == DAE_ALL_NONNEGATIVE)
    {
      for (int i = 0; i < num_local_values; ++i)
        NV_Ith(array, i) = 1.0;
    }
    else if (constraints == DAE_ALL_POSITIVE)
    {
      for (int i = 0; i < num_local_values; ++i)
        NV_Ith(array, i) = 2.0;
    }
    else
    {
      for (int i = 0; i < num_local_values; ++i)
      {
        real_t constraint = constraints[i];
        if (constraint == DAE_UNCONSTRAINED)
          NV_Ith(array, i) = 0.0;
        else if (constraint == DAE_NEGATIVE)
          NV_Ith(array, i) = -2.0;
        else if (constraint == DAE_NONPOSITIVE)
          NV_Ith(array, i) = -1.0;
        else if (constraint == DAE_NONNEGATIVE)
          NV_Ith(array, i) = 1.0;
        else // (constraint == DAE_POSITIVE)
          NV_Ith(array, i) = 2.0;
      }
    }
    IDASetConstraints(integ->ida, array);
  }
  N_VDestroy(array);

  // If we have algebraic equations, we suppress their contributions to 
  // the estimate of the truncation error.
  if (has_algebraic_eqns)
    IDASetSuppressAlg(integ->ida, 1);

  // We use modified Gram-Schmidt orthogonalization.
  IDASpilsSetGSType(integ->ida, MODIFIED_GS);

  // Set up the Jacobian function if given.
  if (integ->Jy != NULL)
    IDASpilsSetJacTimesVecFn(integ->ida, eval_Jy);

  // Set up preconditioner machinery.
  IDASpilsSetPreconditioner(integ->ida, set_up_preconditioner,
                            solve_preconditioner_system);

  // Set some default tolerances:
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  dae_integrator_set_tolerances(integ, 1e-4, 1.0);

  // Set up a maximum number of steps to take during the integration.
//  IDASetMaxNumSteps(integ->ida, 500); // default is 500.

  return integ;
}
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

void dae_integrator_free(dae_integrator_t* integ)
{
  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the IDA stuff.
  N_VDestroy(integ->U_dot);
  polymec_free(integ->U_with_ghosts);
  N_VDestroy(integ->U);
  polymec_free(integ->U_dot_with_ghosts);
  IDAFree(&integ->ida);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  if (integ->error_weights != NULL)
    polymec_free(integ->error_weights);
  polymec_free(integ);
}

void* dae_integrator_context(dae_integrator_t* integ)
{
  return integ->context;
}

int dae_integrator_order(dae_integrator_t* integ)
{
  return integ->order;
}

newton_pc_t* dae_integrator_preconditioner(dae_integrator_t* integrator)
{
  return integrator->precond;
}

void dae_integrator_set_tolerances(dae_integrator_t* integ,
                                   real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  // Clear any existing error weight function.
  integ->compute_weights = NULL;
  if (integ->error_weights != NULL)
  {
    polymec_free(integ->error_weights);
    integ->error_weights = NULL;
  }


  // Set the tolerances.
  IDASStolerances(integ->ida, relative_tol, absolute_tol);
}

// Constant error weight adaptor function.
static void use_constant_weights(void* context, real_t* y, real_t* weights)
{
  dae_integrator_t* integ = context;
  ASSERT(integ->error_weights != NULL);
  memcpy(weights, integ->error_weights, sizeof(real_t) * integ->num_local_values);
}

void dae_integrator_set_error_weights(dae_integrator_t* integ, real_t* weights)
{
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
  dae_integrator_set_error_weight_function(integ, use_constant_weights);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  dae_integrator_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));
  return 0;
}

void dae_integrator_set_error_weight_function(dae_integrator_t* integ,
                                              void (*compute_weights)(void* context, real_t* y, real_t* weights))
{
  ASSERT(compute_weights != NULL);
  integ->compute_weights = compute_weights;
  IDAWFtolerances(integ->ida, compute_error_weights);
}

void dae_integrator_eval_residual(dae_integrator_t* integ, real_t t, real_t* X, real_t* X_dot, real_t* F);
void dae_integrator_eval_residual(dae_integrator_t* integ, real_t t, real_t* X, real_t* X_dot, real_t* F)
{
  START_FUNCTION_TIMER();
  memcpy(integ->U_with_ghosts, X, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, X_dot, sizeof(real_t) * integ->num_local_values);
  integ->F(integ->context, t, integ->U_with_ghosts, integ->U_dot_with_ghosts, F);
  STOP_FUNCTION_TIMER();
}

void dae_integrator_set_max_dt(dae_integrator_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0);
  integ->max_dt = max_dt;
  IDASetMaxStep(integ->ida, max_dt);
}

void dae_integrator_set_stop_time(dae_integrator_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
  IDASetStopTime(integ->ida, stop_time);
}

bool dae_integrator_step(dae_integrator_t* integ, real_t max_dt, real_t* t, real_t* X, real_t* X_dot)
{
  START_FUNCTION_TIMER();
  int status = IDA_SUCCESS;

  // if we haven't been initialized, we need to copy in the data.
  if (!integ->initialized)
    dae_integrator_reset(integ, *t, X, X_dot, false);

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > integ->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = IDASolve(integ->ida, t2, &integ->t, integ->U, integ->U_dot, IDA_ONE_STEP);
    if ((status != IDA_SUCCESS) && (status != IDA_TSTOP_RETURN))
    {
      integ->status_message = get_status_message(status, integ->t);
      STOP_FUNCTION_TIMER();
      return false;
    }

    if ((t2 - *t) < (integ->t - *t))
      log_detail("dae_integrator: took internal step dt = %g", integ->t - *t);
  }

  // If we integrated past t2, interpolate to t2.
  if (integ->t > t2)
  {
    status = IDAGetDky(integ->ida, t2, 0, integ->U);
    status = IDAGetDky(integ->ida, t2, 1, integ->U_dot);
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
  if ((status == IDA_SUCCESS) || (status == IDA_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(X, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    memcpy(X_dot, NV_DATA(integ->U_dot), sizeof(real_t) * integ->num_local_values); 
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

void dae_integrator_reset(dae_integrator_t* integ, 
                          real_t t, 
                          real_t* X, 
                          real_t* X_dot,
                          bool correct_initial_conditions)
{
  // Reset the preconditioner.
  newton_pc_reset(integ->precond, t);

  // Reset the integrator itself.
  integ->t = t;
  memcpy(NV_DATA(integ->U), X, sizeof(real_t) * integ->num_local_values); 
  memcpy(NV_DATA(integ->U_dot), X_dot, sizeof(real_t) * integ->num_local_values); 
  IDAReInit(integ->ida, integ->t, integ->U, integ->U_dot);
  integ->initialized = true;

  if (correct_initial_conditions)
  {
    // Correct the initial conditions for X given those for X_dot.
    log_detail("dae_integrator: correcting initial conditions...");
    int err = IDACalcIC(integ->ida, IDA_Y_INIT, integ->max_dt);
    if (err != IDA_SUCCESS)
    {
      log_detail("dae_integrator: could not correct initial conditions:");
      if (err == IDA_LSETUP_FAIL)
        log_detail("  linear solver setup failed unrecoverably.");
      else if (err == IDA_LINIT_FAIL)
        log_detail("  linear solver initialization failed.");
      else if (err == IDA_LSOLVE_FAIL)
        log_detail("  linear solve failed.");
      else if (err == IDA_BAD_EWT)
        log_detail("  zero component of error weight function due to IC correction.");
      else if (err == IDA_FIRST_RES_FAIL)
        log_detail("  first call to residual failed, causing IC correction to fail.");
      else if (err == IDA_RES_FAIL)
        log_detail("  call to residual failed unrecoverably.");
      else if (err == IDA_CONSTR_FAIL)
        log_detail("  IC corrections were unable to meet constraints.");
      else if (err == IDA_LINESEARCH_FAIL)
        log_detail("  IC correction failed in line search.");
      else
      {
        log_detail("  IC correction failed to converge.");
      }
    }
    IDAGetConsistentIC(integ->ida, integ->U, integ->U_dot);
    memcpy(X, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    memcpy(X_dot, NV_DATA(integ->U_dot), sizeof(real_t) * integ->num_local_values); 
  }

  // Write out debugging info.
  if (log_level() == LOG_DEBUG)
  {
    long num_backtracks;
    IDAGetNumBacktrackOps(integ->ida, &num_backtracks);
    log_debug("dae_integrator: backtracked %ld times in IC correction.", num_backtracks);
  }
}

void dae_integrator_get_diagnostics(dae_integrator_t* integ, 
                                    dae_integrator_diagnostics_t* diagnostics)
{
  diagnostics->status_message = integ->status_message; // borrowed!
  IDAGetNumSteps(integ->ida, &diagnostics->num_steps);
  IDAGetLastOrder(integ->ida, &diagnostics->order_of_last_step);
  IDAGetActualInitStep(integ->ida, &diagnostics->initial_step_size);
  IDAGetLastStep(integ->ida, &diagnostics->last_step_size);
  IDAGetNumResEvals(integ->ida, &diagnostics->num_residual_evaluations);
  IDAGetNumLinSolvSetups(integ->ida, &diagnostics->num_linear_solve_setups);
  IDAGetNumErrTestFails(integ->ida, &diagnostics->num_error_test_failures);
  IDAGetNumNonlinSolvIters(integ->ida, &diagnostics->num_nonlinear_solve_iterations);
  IDAGetNumNonlinSolvConvFails(integ->ida, &diagnostics->num_nonlinear_solve_convergence_failures);
  IDASpilsGetNumLinIters(integ->ida, &diagnostics->num_linear_solve_iterations);
  IDASpilsGetNumPrecEvals(integ->ida, &diagnostics->num_preconditioner_evaluations);
  IDASpilsGetNumPrecSolves(integ->ida, &diagnostics->num_preconditioner_solves);
  IDASpilsGetNumConvFails(integ->ida, &diagnostics->num_linear_solve_convergence_failures);
}

void dae_integrator_diagnostics_fprintf(dae_integrator_diagnostics_t* diagnostics, 
                                        FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "dae integrator diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Order of last step: %d\n", diagnostics->order_of_last_step);
  fprintf(stream, "  Initial step size: %g\n", diagnostics->initial_step_size);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Num residual evaluations: %d\n", (int)diagnostics->num_residual_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
}


