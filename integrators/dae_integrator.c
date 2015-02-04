// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "integrators/dae_integrator.h"

#include "ida/ida.h"
#include "ida/ida_spils.h"
#include "ida/ida_spgmr.h"
#include "ida/ida_spbcgs.h"
#include "ida/ida_sptfqmr.h"

// Types of linear solvers.
typedef enum 
{
  GMRES,
  BICGSTAB,
  TFQMR
} solver_type_t;

struct dae_integrator_t 
{
  char* name;
  void* context;
  dae_integrator_vtable vtable;
  int order;
  MPI_Comm comm;
  solver_type_t solver_type;
  bool initialized;

  int N; // dimension of system.

  // IDA data structures.
  void* ida;
  N_Vector x, x_dot;
  real_t current_time;
  int max_krylov_dim;
  char* status_message; // status of most recent integration.
  real_t max_dt, stop_time;

  // Error weight function.
  dae_integrator_error_weight_func compute_weights;

  // Preconditioning stuff. 
  newton_pc_t* precond;
};

// This function wraps around the user-supplied right hand side.
static int evaluate_residual(real_t t, N_Vector x, N_Vector x_dot, 
                             N_Vector F, void* context)
{
  dae_integrator_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);
  real_t* Fx = NV_DATA(F);

  // Evaluate the residual.
  return integ->vtable.residual(integ->context, t, xx, xxd, Fx);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector x_dot, N_Vector F,
                                 real_t cj, void* context,
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  dae_integrator_t* integ = context;

  // Form the linear combination dFdx + cj * dFdxdot.
  newton_pc_setup(integ->precond, 0.0, 1.0, cj, t, NV_DATA(x), NV_DATA(x_dot));

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
  dae_integrator_t* integ = context;
  
  // FIXME: Apply scaling if needed.
  // Copy the contents of the residual to the output vector.
  memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * integ->N);

  // Solve it.
  if (newton_pc_solve(integ->precond, NV_DATA(z)))
    return 0;
  else 
    return 1; // recoverable error.
}

static dae_integrator_t* dae_integrator_new(const char* name, 
                                            void* context,
                                            MPI_Comm comm,
                                            int N,
                                            dae_integrator_vtable vtable,
                                            int order,
                                            solver_type_t solver_type,
                                            int max_krylov_dim)
{
  ASSERT(N > 0);
  ASSERT(order > 0);
  ASSERT(order <= 5);
  ASSERT(vtable.residual != NULL);
  ASSERT(max_krylov_dim >= 3);

  dae_integrator_t* integ = polymec_malloc(sizeof(dae_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->comm = comm;
  integ->vtable = vtable;
  integ->order = order;
  integ->solver_type = solver_type;
  integ->current_time = 0.0;
  integ->N = N;
  integ->max_krylov_dim = max_krylov_dim;
  integ->initialized = false;
  integ->max_dt = FLT_MAX;
  integ->status_message = NULL;

  // Set up KINSol and accessories.
  integ->x = N_VNew(comm, N);
  integ->x_dot = N_VNew(comm, N);
  integ->ida = IDACreate();
  IDASetMaxOrd(integ->ida, integ->order);
  IDASetUserData(integ->ida, integ);
  IDAInit(integ->ida, evaluate_residual, integ->current_time, integ->x, integ->x_dot);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver_type == GMRES)
    IDASpgmr(integ->ida, max_krylov_dim); 
  else if (solver_type == BICGSTAB)
    IDASpbcg(integ->ida, max_krylov_dim);
  else
    IDASptfqmr(integ->ida, max_krylov_dim);

  // We use modified Gram-Schmidt orthogonalization.
  IDASpilsSetGSType(integ->ida, MODIFIED_GS);

  // Set up preconditioner machinery.
  IDASpilsSetPreconditioner(integ->ida, set_up_preconditioner,
                           solve_preconditioner_system);
  integ->precond = NULL;

  // Algebraic constraints.
  if (integ->vtable.set_constraints != NULL)
  {
    N_Vector constraints = N_VNew(integ->comm, N);
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

dae_integrator_t* gmres_dae_integrator_new(const char* name,
                                           void* context,
                                           MPI_Comm comm,
                                           int N,
                                           dae_integrator_vtable vtable,
                                           int order,
                                           int max_krylov_dim)
{
  return dae_integrator_new(name, context, comm, N, vtable, order, GMRES, 
                            max_krylov_dim);
}

dae_integrator_t* bicgstab_dae_integrator_new(const char* name,
                                              void* context,
                                              MPI_Comm comm,
                                              int N,
                                              dae_integrator_vtable vtable,
                                              int order,
                                              int max_krylov_dim)
{
  return dae_integrator_new(name, context, comm, N, vtable, order, BICGSTAB, 
                            max_krylov_dim);
}

dae_integrator_t* tfqmr_dae_integrator_new(const char* name,
                                           void* context,
                                           MPI_Comm comm,
                                           int N,
                                           dae_integrator_vtable vtable,
                                           int order,
                                           int max_krylov_dim)
{
  return dae_integrator_new(name, context, comm, N, vtable, order, TFQMR, 
                            max_krylov_dim);
}

void dae_integrator_free(dae_integrator_t* integ)
{
  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the IDA stuff.
  N_VDestroy(integ->x_dot);
  N_VDestroy(integ->x);
  IDAFree(&integ->ida);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  polymec_free(integ->name);
  polymec_free(integ);
}

char* dae_integrator_name(dae_integrator_t* integ)
{
  return integ->name;
}

void* dae_integrator_context(dae_integrator_t* integ)
{
  return integ->context;
}

int dae_integrator_order(dae_integrator_t* integ)
{
  return integ->order;
}

void dae_integrator_set_preconditioner(dae_integrator_t* integrator,
                                       newton_pc_t* precond)
{
  integrator->precond = precond;
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
  integ->vtable.residual(integ->context, t, X, X_dot, F);
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

bool dae_integrator_step(dae_integrator_t* integ, real_t* t, real_t* X, real_t* X_dot)
{
  // Copy in the solution and its time derivative.
  memcpy(NV_DATA(integ->x), X, sizeof(real_t) * integ->N); 
  memcpy(NV_DATA(integ->x_dot), X_dot, sizeof(real_t) * integ->N); 

  if (!integ->initialized)
  {
    integ->current_time = *t;
    IDAReInit(integ->ida, integ->current_time, integ->x, integ->x_dot);
    integ->initialized = true;
  }
  else if (fabs(integ->current_time - *t) > 1e-14)
  {
    // Reset the integrator if t != current_time.
    integ->current_time = *t;
    IDAReInit(integ->ida, integ->current_time, integ->x, integ->x_dot);
  }

  // Integrate.
  real_t t2 = *t + integ->max_dt;
  int status = IDASolve(integ->ida, t2, &integ->current_time, integ->x, integ->x_dot, IDA_ONE_STEP);
  
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
    *t = integ->current_time;
    memcpy(X, NV_DATA(integ->x), sizeof(real_t) * integ->N); 
    memcpy(X_dot, NV_DATA(integ->x_dot), sizeof(real_t) * integ->N); 
    return true;
  }
  else
  {
    if (status == IDA_TOO_MUCH_WORK)
    {
      char err[1024];
      snprintf(err, 1024, "Integrator stopped at t = %g after maximum number of steps.", integ->current_time);
      integ->status_message = string_dup(err);
    }
    else if (status == IDA_TOO_MUCH_ACC)
      integ->status_message = string_dup("Integrator could not achieve desired level of accuracy.");
    else if (status == IDA_ERR_FAIL)
      integ->status_message = string_dup("Integrator encountered too many error test failures.");
    else if (status == IDA_CONV_FAIL)
      integ->status_message = string_dup("Integrator encountered too many convergence test failures.");
    else if (status == IDA_LINIT_FAIL)
      integ->status_message = string_dup("Integrator's linear solver failed to initialize.");
    else if (status == IDA_LSETUP_FAIL)
      integ->status_message = string_dup("Integrator's linear solver setup failed.");
    else if (status == IDA_LSOLVE_FAIL)
      integ->status_message = string_dup("Integrator's linear solver failed.");
    else if (status == IDA_CONSTR_FAIL)
      integ->status_message = string_dup("Integrator's inequality constraints were violated unrecoverably.");
    else if (status == IDA_RES_FAIL)
      integ->status_message = string_dup("Integrator's residual function failed unrecoverably.");
    else if (status == IDA_REP_RES_ERR)
      integ->status_message = string_dup("Integrator encountered too many recoverable RHS failures.");
    else if (status == IDA_RTFUNC_FAIL)
      integ->status_message = string_dup("Integrator encountered a failure in the rootfinding function.");
    return false;
  }
}

void dae_integrator_reset(dae_integrator_t* integ, real_t t, real_t* X, real_t* X_dot)
{
  integ->current_time = t;
  memcpy(NV_DATA(integ->x), X, sizeof(real_t) * integ->N); 
  memcpy(NV_DATA(integ->x_dot), X_dot, sizeof(real_t) * integ->N); 
  IDAReInit(integ->ida, integ->current_time, integ->x, integ->x_dot);
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


