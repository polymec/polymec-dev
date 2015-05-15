// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

struct dae_integrator_t 
{
  char* name;
  void* context;
  dae_integrator_vtable vtable;
  int order;
  MPI_Comm comm;
  dae_krylov_t solver_type;
  bool initialized;

  int num_local_values, num_remote_values;

  // IDA data structures.
  void* ida;
  N_Vector x, x_dot;
  real_t* x_with_ghosts;
  real_t* x_dot_with_ghosts;
  int max_krylov_dim;
  char* status_message; // status of most recent integration.
  real_t max_dt, stop_time;

  // Current simulation time;
  real_t t;

  // Error weight function.
  dae_integrator_error_weight_func compute_weights;

  // Preconditioning stuff. 
  newton_pc_t* precond;
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
static int evaluate_residual(real_t t, N_Vector x, N_Vector x_dot, 
                             N_Vector F, void* context)
{
  dae_integrator_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);
  real_t* Fx = NV_DATA(F);

  // Evaluate the residual using vectors with ghosts.
  memcpy(integ->x_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->x_dot_with_ghosts, xxd, sizeof(real_t) * integ->num_local_values);
  return integ->vtable.residual(integ->context, t, integ->x_with_ghosts, 
                                integ->x_dot_with_ghosts, Fx);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector x_dot, N_Vector F,
                                 real_t cj, void* context,
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  START_FUNCTION_TIMER();
  dae_integrator_t* integ = context;

  // Form the linear combination dFdx + cj * dFdxdot.
  memcpy(integ->x_with_ghosts, NV_DATA(x), sizeof(real_t) * integ->num_local_values);
  memcpy(integ->x_dot_with_ghosts, NV_DATA(x_dot), sizeof(real_t) * integ->num_local_values);
  newton_pc_setup(integ->precond, 0.0, 1.0, cj, t, integ->x_with_ghosts, integ->x_dot_with_ghosts);
  STOP_FUNCTION_TIMER();
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector x, N_Vector x_dot, 
                                       N_Vector F, N_Vector r, N_Vector z, 
                                       real_t cj, real_t delta, 
                                       void* context, N_Vector work)
{
  START_FUNCTION_TIMER();
  dae_integrator_t* integ = context;
  
  // FIXME: Apply scaling if needed.
  // Copy the contents of the residual to the output vector.
  memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * integ->num_local_values);

  // Solve it.
  int result = (newton_pc_solve(integ->precond, NV_DATA(z))) ? 0 : 1;
  STOP_FUNCTION_TIMER();
  return result;
}

dae_integrator_t* dae_integrator_new(int order,
                                     MPI_Comm comm,
                                     dae_equation_t* equation_types,
                                     int num_local_values,
                                     int num_remote_values,
                                     void* context,
                                     dae_integrator_vtable vtable,
                                     newton_pc_t* precond,
                                     dae_krylov_t solver_type,
                                     int max_krylov_dim)
{
  ASSERT(order > 0);
  ASSERT(order <= 5);
  ASSERT(equation_types != NULL);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(vtable.residual != NULL);
  ASSERT(precond != NULL);
  ASSERT(max_krylov_dim >= 3);

  dae_integrator_t* integ = polymec_malloc(sizeof(dae_integrator_t));
  integ->context = context;
  integ->comm = comm;
  integ->vtable = vtable;
  integ->order = order;
  integ->solver_type = solver_type;
  integ->t = 0.0;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->precond = precond;
  integ->max_krylov_dim = max_krylov_dim;
  integ->initialized = false;
  integ->max_dt = FLT_MAX;
  integ->status_message = NULL;

  // Set up IDA and accessories.
  integ->x = N_VNew(comm, num_local_values);
  integ->x_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->x_dot = N_VNew(comm, num_local_values);
  integ->x_dot_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  memset(NV_DATA(integ->x), 0, sizeof(real_t) * num_local_values);
  memset(NV_DATA(integ->x_dot), 0, sizeof(real_t) * num_local_values);

  integ->ida = IDACreate();
  IDASetMaxOrd(integ->ida, integ->order);
  IDASetUserData(integ->ida, integ);
  IDAInit(integ->ida, evaluate_residual, integ->t, integ->x, integ->x_dot);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver_type == DAE_GMRES)
    IDASpgmr(integ->ida, max_krylov_dim); 
  else if (solver_type == DAE_BICGSTAB)
    IDASpbcg(integ->ida, max_krylov_dim);
  else
    IDASptfqmr(integ->ida, max_krylov_dim);

  // We set the equation types.
  N_Vector types = N_VNew(comm, num_local_values);
  bool has_algebraic_eqns = false;
  for (int i = 0; i < num_local_values; ++i)
  {
    if (equation_types[i] == DAE_ALGEBRAIC)
    {
      has_algebraic_eqns = true;
      NV_Ith(types, i) = 0.0;
    }
    else
      NV_Ith(types, i) = 1.0;
  }
  IDASetId(integ->ida, types);
  N_VDestroy(types);

  // If we have algebraic equations, we suppress their contributions to 
  // the estimate of the truncation error.
  if (has_algebraic_eqns)
    IDASetSuppressAlg(integ->ida, 1);

  // We use modified Gram-Schmidt orthogonalization.
  IDASpilsSetGSType(integ->ida, MODIFIED_GS);

  // Set up preconditioner machinery.
  IDASpilsSetPreconditioner(integ->ida, set_up_preconditioner,
                           solve_preconditioner_system);

  // Algebraic constraints.
  if (integ->vtable.set_constraints != NULL)
  {
    N_Vector constraints = N_VNew(integ->comm, num_local_values);
    integ->vtable.set_constraints(integ->context, NV_DATA(constraints));
    IDASetConstraints(integ->ida, constraints);
    N_VDestroy(constraints);
  }

  // Set some default tolerances:
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  dae_integrator_set_tolerances(integ, 1e-4, 1.0);

  // Set up a maximum number of steps to take during the integration.
//  IDASetMaxNumSteps(integ->ida, 500); // default is 500.

  return integ;
}

void dae_integrator_free(dae_integrator_t* integ)
{
  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the IDA stuff.
  N_VDestroy(integ->x_dot);
  polymec_free(integ->x_with_ghosts);
  N_VDestroy(integ->x);
  polymec_free(integ->x_dot_with_ghosts);
  IDAFree(&integ->ida);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
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

void dae_integrator_set_tolerances(dae_integrator_t* integrator,
                                   real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  // Clear any existing error weight function.
  integrator->compute_weights = NULL;

  // Set the tolerances.
  IDASStolerances(integrator->ida, relative_tol, absolute_tol);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  dae_integrator_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));
  return 0;
}

void dae_integrator_set_error_weight_function(dae_integrator_t* integrator,
                                              dae_integrator_error_weight_func compute_weights)
{
  ASSERT(compute_weights != NULL);
  integrator->compute_weights = compute_weights;
  IDAWFtolerances(integrator->ida, compute_error_weights);
}

void dae_integrator_eval_residual(dae_integrator_t* integ, real_t t, real_t* X, real_t* X_dot, real_t* F)
{
  START_FUNCTION_TIMER();
  memcpy(integ->x_with_ghosts, X, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->x_dot_with_ghosts, X_dot, sizeof(real_t) * integ->num_local_values);
  integ->vtable.residual(integ->context, t, integ->x_with_ghosts, integ->x_dot_with_ghosts, F);
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
    dae_integrator_reset(integ, *t, X, X_dot);

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > integ->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = IDASolve(integ->ida, t2, &integ->t, integ->x, integ->x_dot, IDA_ONE_STEP);
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
    status = IDAGetDky(integ->ida, t2, 0, integ->x);
    status = IDAGetDky(integ->ida, t2, 1, integ->x_dot);
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
    memcpy(X, NV_DATA(integ->x), sizeof(real_t) * integ->num_local_values); 
    memcpy(X_dot, NV_DATA(integ->x_dot), sizeof(real_t) * integ->num_local_values); 
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

void dae_integrator_reset(dae_integrator_t* integ, real_t t, real_t* X, real_t* X_dot)
{
  integ->t = t;
  memcpy(NV_DATA(integ->x), X, sizeof(real_t) * integ->num_local_values); 
  memcpy(NV_DATA(integ->x_dot), X_dot, sizeof(real_t) * integ->num_local_values); 
  IDAReInit(integ->ida, integ->t, integ->x, integ->x_dot);
  integ->initialized = true;
}

void dae_integrator_get_diagnostics(dae_integrator_t* integrator, 
                                    dae_integrator_diagnostics_t* diagnostics)
{
  diagnostics->status_message = integrator->status_message; // borrowed!
  IDAGetNumSteps(integrator->ida, &diagnostics->num_steps);
  IDAGetLastOrder(integrator->ida, &diagnostics->order_of_last_step);
  IDAGetActualInitStep(integrator->ida, &diagnostics->initial_step_size);
  IDAGetLastStep(integrator->ida, &diagnostics->last_step_size);
  IDAGetNumResEvals(integrator->ida, &diagnostics->num_residual_evaluations);
  IDAGetNumLinSolvSetups(integrator->ida, &diagnostics->num_linear_solve_setups);
  IDAGetNumErrTestFails(integrator->ida, &diagnostics->num_error_test_failures);
  IDAGetNumNonlinSolvIters(integrator->ida, &diagnostics->num_nonlinear_solve_iterations);
  IDAGetNumNonlinSolvConvFails(integrator->ida, &diagnostics->num_nonlinear_solve_convergence_failures);
  IDASpilsGetNumLinIters(integrator->ida, &diagnostics->num_linear_solve_iterations);
  IDASpilsGetNumPrecEvals(integrator->ida, &diagnostics->num_preconditioner_evaluations);
  IDASpilsGetNumPrecSolves(integrator->ida, &diagnostics->num_preconditioner_solves);
  IDASpilsGetNumConvFails(integrator->ida, &diagnostics->num_linear_solve_convergence_failures);
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


