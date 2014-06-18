// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/sundials_helpers.h"
#include "integrators/ode_integrator.h"

#include "cvode/cvode.h"
#include "cvode/cvode_spils.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"
#include "cvode/cvode_sptfqmr.h"

// Types of linear solvers.
typedef enum 
{
  GMRES,
  BICGSTAB,
  TFQMR
} solver_type_t;

struct ode_integrator_t 
{
  char* name;
  void* context;
  ode_integrator_vtable vtable;
  int order;
  MPI_Comm comm;
  solver_type_t solver_type;
  bool initialized;

  int N; // dimension of system.

  // CVODE data structures.
  void* cvode;
  N_Vector x; 
  real_t current_time;
  int max_krylov_dim;
  char* status_message; // status of most recent integration.
  real_t max_dt, stop_time;

  // Error weight function.
  ode_integrator_error_weight_func compute_weights;

  // Preconditioning stuff.
  preconditioner_t* precond;
};

// This function wraps around the user-supplied right hand side.
static int evaluate_rhs(real_t t, N_Vector x, N_Vector x_dot, void* context)
{
  ode_integrator_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);

  // Evaluate the RHS.
  return integ->vtable.rhs(integ->context, t, xx, xxd);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  ode_integrator_t* integ = context;
  if (!jacobian_is_current)
  {
    preconditioner_setup(integ->precond, 1.0, -gamma, 0.0, t, NV_DATA(x), NULL);
    *jacobian_was_updated = 1;
  }
  else
    *jacobian_was_updated = 0;
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector x, N_Vector F, 
                                       N_Vector r, N_Vector z, 
                                       real_t gamma, real_t delta, 
                                       int lr, void* context, 
                                       N_Vector work)
{
  ASSERT(lr == 1); // Left preconditioning only.

  ode_integrator_t* integ = context;
  
  // FIXME: Apply scaling if needed.
  // Copy the contents of the RHS to the output vector.
  memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * integ->N);

  // Solve it.
  if (preconditioner_solve(integ->precond, NV_DATA(z)))
    return 0;
  else 
    return 1; // recoverable error.
}

static ode_integrator_t* ode_integrator_new(const char* name, 
                                            void* context,
                                            MPI_Comm comm,
                                            int N,
                                            ode_integrator_vtable vtable,
                                            int order,
                                            solver_type_t solver_type,
                                            int max_krylov_dim)
{
  ASSERT(N > 0);
  ASSERT(order > 0);
  ASSERT(order <= 5);
  ASSERT(vtable.rhs != NULL);
  ASSERT(max_krylov_dim >= 3);

  ode_integrator_t* integ = polymec_malloc(sizeof(ode_integrator_t));
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
  integ->stop_time = FLT_MAX;
  integ->status_message = NULL;

  // Set up KINSol and accessories.
  integ->x = N_VNew(comm, N);
  integ->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetMaxOrd(integ->cvode, integ->order);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, evaluate_rhs, integ->current_time, integ->x);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver_type == GMRES)
  {
    CVSpgmr(integ->cvode, PREC_LEFT, max_krylov_dim); 
    // We use modified Gram-Schmidt orthogonalization.
    CVSpilsSetGSType(integ->cvode, MODIFIED_GS);
  }
  else if (solver_type == BICGSTAB)
    CVSpbcg(integ->cvode, PREC_LEFT, max_krylov_dim);
  else
    CVSptfqmr(integ->cvode, PREC_LEFT, max_krylov_dim);

  CVSpilsSetPreconditioner(integ->cvode, set_up_preconditioner,
                           solve_preconditioner_system);

  integ->precond = NULL;

  // Set some default tolerances:
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  ode_integrator_set_tolerances(integ, 1e-4, 1.0);

  // Set up a maximum number of steps to take during the integration.
//  CVodeSetMaxNumSteps(integ->cvode, 500); // default is 500.

  return integ;
}

ode_integrator_t* gmres_ode_integrator_new(const char* name,
                                           void* context,
                                           MPI_Comm comm,
                                           int N,
                                           ode_integrator_vtable vtable,
                                           int order,
                                           int max_krylov_dim)
{
  return ode_integrator_new(name, context, comm, N, vtable, order, GMRES, 
                            max_krylov_dim);
}

ode_integrator_t* bicgstab_ode_integrator_new(const char* name,
                                              void* context,
                                              MPI_Comm comm,
                                              int N,
                                              ode_integrator_vtable vtable,
                                              int order,
                                              int max_krylov_dim)
{
  return ode_integrator_new(name, context, comm, N, vtable, order, BICGSTAB, 
                            max_krylov_dim);
}

ode_integrator_t* tfqmr_ode_integrator_new(const char* name,
                                           void* context,
                                           MPI_Comm comm,
                                           int N,
                                           ode_integrator_vtable vtable,
                                           int order,
                                           int max_krylov_dim)
{
  return ode_integrator_new(name, context, comm, N, vtable, order, TFQMR, 
                            max_krylov_dim);
}

void ode_integrator_free(ode_integrator_t* integ)
{
  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    preconditioner_free(integ->precond);

  // Kill the CVode stuff.
  N_VDestroy(integ->x);
  CVodeFree(&integ->cvode);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  polymec_free(integ->name);
  polymec_free(integ);
}

char* ode_integrator_name(ode_integrator_t* integ)
{
  return integ->name;
}

void* ode_integrator_context(ode_integrator_t* integ)
{
  return integ->context;
}

int ode_integrator_order(ode_integrator_t* integ)
{
  return integ->order;
}

void ode_integrator_set_preconditioner(ode_integrator_t* integrator,
                                       preconditioner_t* precond)
{
  integrator->precond = precond;
}

preconditioner_t* ode_integrator_preconditioner(ode_integrator_t* integrator)
{
  return integrator->precond;
}

void ode_integrator_set_stability_limit_detection(ode_integrator_t* integrator,
                                                   bool use_detection)
{
  CVodeSetStabLimDet(integrator->cvode, use_detection);
}

void ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                    real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  // Clear any existing error weight function.
  integrator->compute_weights = NULL;

  // Set the tolerances.
  CVodeSStolerances(integrator->cvode, relative_tol, absolute_tol);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  ode_integrator_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));
  return 0;
}

void ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                              ode_integrator_error_weight_func compute_weights)
{
  ASSERT(compute_weights != NULL);
  integrator->compute_weights = compute_weights;
  CVodeWFtolerances(integrator->cvode, compute_error_weights);
}

void ode_integrator_eval_rhs(ode_integrator_t* integ, real_t t, real_t* X, real_t* rhs)
{
  integ->vtable.rhs(integ->context, t, X, rhs);
}

void ode_integrator_set_max_dt(ode_integrator_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0);
  integ->max_dt = max_dt;
  CVodeSetMaxStep(integ->cvode, max_dt);
}

void ode_integrator_set_stop_time(ode_integrator_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
  CVodeSetStopTime(integ->cvode, stop_time);
}

bool ode_integrator_step(ode_integrator_t* integ, real_t* t, real_t* X)
{
  // Copy in the solution.
  memcpy(NV_DATA(integ->x), X, sizeof(real_t) * integ->N); 

  if (!integ->initialized)
  {
    integ->current_time = *t;
    CVodeReInit(integ->cvode, integ->current_time, integ->x);
    integ->initialized = true;
  }
  else if (fabs(integ->current_time - *t) > 1e-14)
  {
    // Reset the integrator if t != current_time.
    integ->current_time = *t;
    CVodeReInit(integ->cvode, integ->current_time, integ->x);
  }

  // Integrate.
  real_t t2 = *t + integ->max_dt;
  int status = CVode(integ->cvode, t2, integ->x, &integ->current_time, CV_ONE_STEP);
  
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
    *t = integ->current_time;
    memcpy(X, NV_DATA(integ->x), sizeof(real_t) * integ->N); 
    return true;
  }
  else
  {
    if (status == CV_TOO_CLOSE)
      integ->status_message = string_dup("t1 and t2 are too close to each other.");
    else if (status == CV_TOO_MUCH_WORK)
    {
      char err[1024];
      snprintf(err, 1024, "Integrator stopped at t = %g after maximum number of steps.", integ->current_time);
      integ->status_message = string_dup(err);
    }
    else if (status == CV_TOO_MUCH_ACC)
      integ->status_message = string_dup("Integrator could not achieve desired level of accuracy.");
    else if (status == CV_ERR_FAILURE)
      integ->status_message = string_dup("Integrator encountered too many error test failures.");
    else if (status == CV_CONV_FAILURE)
      integ->status_message = string_dup("Integrator encountered too many convergence test failures.");
    else if (status == CV_LINIT_FAIL)
      integ->status_message = string_dup("Integrator's linear solver failed to initialize.");
    else if (status == CV_LSETUP_FAIL)
      integ->status_message = string_dup("Integrator's linear solver setup failed.");
    else if (status == CV_LSOLVE_FAIL)
      integ->status_message = string_dup("Integrator's linear solver failed.");
    else if (status == CV_RHSFUNC_FAIL)
      integ->status_message = string_dup("Integrator's RHS function failed unrecoverably.");
    else if (status == CV_FIRST_RHSFUNC_ERR)
      integ->status_message = string_dup("Integrator's first call to RHS function failed.");
    else if (status == CV_REPTD_RHSFUNC_ERR)
      integ->status_message = string_dup("Integrator encountered too many recoverable RHS failures.");
    else if (status == CV_UNREC_RHSFUNC_ERR)
      integ->status_message = string_dup("Integrator failed to recover from a recoverable RHS failure.");
    else if (status == CV_RTFUNC_FAIL)
      integ->status_message = string_dup("Integrator encountered a failure in the rootfinding function.");
    return false;
  }
}

void ode_integrator_reset(ode_integrator_t* integ, real_t t, real_t* X)
{
  integ->current_time = t;
  memcpy(NV_DATA(integ->x), X, sizeof(real_t) * integ->N); 
  CVodeReInit(integ->cvode, integ->current_time, integ->x);
  integ->initialized = true;
}

void ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                    ode_integrator_diagnostics_t* diagnostics)
{
  diagnostics->status_message = integrator->status_message; // borrowed!
  CVodeGetNumSteps(integrator->cvode, &diagnostics->num_steps);
  CVodeGetLastOrder(integrator->cvode, &diagnostics->order_of_last_step);
  CVodeGetLastStep(integrator->cvode, &diagnostics->last_step_size);
  CVodeGetNumRhsEvals(integrator->cvode, &diagnostics->num_rhs_evaluations);
  CVodeGetNumLinSolvSetups(integrator->cvode, &diagnostics->num_linear_solve_setups);
  CVodeGetNumErrTestFails(integrator->cvode, &diagnostics->num_error_test_failures);
  CVodeGetNumNonlinSolvIters(integrator->cvode, &diagnostics->num_nonlinear_solve_iterations);
  CVodeGetNumNonlinSolvConvFails(integrator->cvode, &diagnostics->num_nonlinear_solve_convergence_failures);
  CVSpilsGetNumLinIters(integrator->cvode, &diagnostics->num_linear_solve_iterations);
  CVSpilsGetNumPrecEvals(integrator->cvode, &diagnostics->num_preconditioner_evaluations);
  CVSpilsGetNumPrecSolves(integrator->cvode, &diagnostics->num_preconditioner_solves);
  CVSpilsGetNumConvFails(integrator->cvode, &diagnostics->num_linear_solve_convergence_failures);
}

void ode_integrator_diagnostics_fprintf(ode_integrator_diagnostics_t* diagnostics, 
                                        FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "ODE integrator diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Order of last step: %d\n", diagnostics->order_of_last_step);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Num RHS evaluations: %d\n", (int)diagnostics->num_rhs_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
}


