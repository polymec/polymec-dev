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
#include "integrators/am_ode_integrator.h"
#include "integrators/nonlinear_preconditioner.h"
#include "cvode/cvode.h"

// JFNK stuff.
#include "cvode/cvode_spils.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"
#include "cvode/cvode_sptfqmr.h"

typedef struct
{
  MPI_Comm comm;
  int N;
  void* context; 
  int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot);
  void (*dtor)(void* context);

  // CVODE data structures.
  bool initialized;
  void* cvode;
  N_Vector x; 
  real_t t;
  char* status_message; // status of most recent integration.

  // JFNK stuff.
  int max_krylov_dim;
  preconditioner_t* precond;
  int (*Jy)(void* context, real_t t, real_t* x, real_t* rhstx, real_t* y, real_t* temp, real_t* Jy);

  // Error weight function.
  am_ode_integrator_error_weight_func compute_weights;
} am_ode_t;

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
static int am_evaluate_rhs(real_t t, N_Vector x, N_Vector x_dot, void* context)
{
  am_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);

  // Evaluate the RHS.
  return integ->rhs(integ->context, t, xx, xxd);
}

static bool am_step(void* context, real_t max_dt, real_t max_t, real_t* t, real_t* x)
{
  am_ode_t* integ = context;
  CVodeSetMaxStep(integ->cvode, max_dt);
  CVodeSetStopTime(integ->cvode, max_t);

  if (!integ->initialized || (fabs(integ->t - *t) > 1e-14))
  {
    CVodeReInit(integ->cvode, *t, integ->x);
    integ->initialized = true;
    integ->t = *t;
  }
  
  // Copy in the solution.
  memcpy(NV_DATA(integ->x), x, sizeof(real_t) * integ->N); 

  // Integrate.
  integ->t = *t;
  real_t t2 = *t + max_dt;
  int status = CVode(integ->cvode, t2, integ->x, &integ->t, CV_ONE_STEP);
  
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
    *t = integ->t;
    memcpy(x, NV_DATA(integ->x), sizeof(real_t) * integ->N); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static bool am_advance(void* context, real_t max_dt, real_t max_t, real_t* x)
{
  am_ode_t* integ = context;
  CVodeSetMaxStep(integ->cvode, max_dt);
  CVodeSetStopTime(integ->cvode, max_t);
  
  // Copy in the solution.
  memcpy(NV_DATA(integ->x), x, sizeof(real_t) * integ->N); 

  // Integrate.
  int status = CVode(integ->cvode, max_t, integ->x, &integ->t, CV_NORMAL);
  
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
    memcpy(x, NV_DATA(integ->x), sizeof(real_t) * integ->N); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static void am_dtor(void* context)
{
  am_ode_t* integ = context;

  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    preconditioner_free(integ->precond);

  // Kill the CVode stuff.
  N_VDestroy(integ->x);
  CVodeFree(&integ->cvode);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  polymec_free(integ);
}

ode_integrator_t* functional_am_ode_integrator_new(int order, 
                                                   MPI_Comm comm, 
                                                   int N,
                                                   void* context, 
                                                   int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                   void (*dtor)(void* context))
{
  ASSERT(order >= 1);
  ASSERT(order <= 12);
  ASSERT(N > 0);
  ASSERT(rhs != NULL);

  am_ode_t* integ = polymec_malloc(sizeof(am_ode_t));
  integ->comm = comm;
  integ->N = N;
  integ->context = context;
  integ->rhs = rhs;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = 0;
  integ->initialized = false;

  // Set up KINSol and accessories.
  integ->x = N_VNew(integ->comm, integ->N);
  integ->cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  CVodeSetMaxOrd(integ->cvode, order);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, am_evaluate_rhs, 0.0, integ->x);

  ode_integrator_vtable vtable = {.step = am_step, .advance = am_advance, .dtor = am_dtor};
  char name[1024];
  snprintf(name, 1024, "Functional Adams-Moulton (order %d)", order);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  am_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  am_ode_t* integ = context;
  if (!jacobian_is_current)
  {
    nonlinear_preconditioner_setup(integ->precond, 1.0, -gamma, 0.0, t, NV_DATA(x), NULL);
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

  am_ode_t* integ = context;
  
  // FIXME: Apply scaling if needed.
  // Copy the contents of the RHS to the output vector.
  memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * integ->N);

  // Solve it.
  if (preconditioner_solve(integ->precond, NV_DATA(z)))
    return 0;
  else 
    return 1; // recoverable error.
}

// Adaptor for J*y function
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector x, N_Vector rhstx, void* context, N_Vector tmp)
{
  am_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* yy = NV_DATA(y);
  real_t* rhs = NV_DATA(rhstx);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);
  return integ->Jy(integ->context, t, xx, rhs, yy, temp, jjy);
}

ode_integrator_t* jfnk_am_ode_integrator_new(int order,
                                             MPI_Comm comm,
                                             int N, 
                                             void* context, 
                                             int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                             int (*Jy)(void* context, real_t t, real_t* x, real_t* rhstx, real_t* y, real_t* temp, real_t* Jy),
                                             void (*dtor)(void* context),
                                             preconditioner_t* precond,
                                             jfnk_am_krylov_t solver_type,
                                             int max_krylov_dim)
{
  ASSERT(order >= 1);
  ASSERT(order <= 12);
  ASSERT(N > 0);
  ASSERT(rhs != NULL);
  ASSERT(precond != NULL);
  ASSERT(max_krylov_dim > 3);

  am_ode_t* integ = polymec_malloc(sizeof(am_ode_t));
  integ->comm = comm;
  integ->N = N;
  integ->context = context;
  integ->rhs = rhs;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = max_krylov_dim;
  integ->initialized = false;

  // Set up KINSol and accessories.
  integ->x = N_VNew(integ->comm, integ->N);
  integ->cvode = CVodeCreate(CV_ADAMS, CV_NEWTON);
  CVodeSetMaxOrd(integ->cvode, order);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, am_evaluate_rhs, 0.0, integ->x);

  // Set up the solver type.
  if (solver_type == JFNK_AM_GMRES)
  {
    CVSpgmr(integ->cvode, PREC_LEFT, max_krylov_dim); 
    // We use modified Gram-Schmidt orthogonalization.
    CVSpilsSetGSType(integ->cvode, MODIFIED_GS);
  }
  else if (solver_type == JFNK_AM_BICGSTAB)
    CVSpbcg(integ->cvode, PREC_LEFT, max_krylov_dim);
  else
    CVSptfqmr(integ->cvode, PREC_LEFT, max_krylov_dim);

  // Set up the Jacobian function and preconditioner.
  if (Jy != NULL)
    CVSpilsSetJacTimesVecFn(integ->cvode, eval_Jy);
  integ->precond = precond;
  CVSpilsSetPreconditioner(integ->cvode, set_up_preconditioner,
                           solve_preconditioner_system);

  ode_integrator_vtable vtable = {.step = am_step, .advance = am_advance, .dtor = am_dtor};
  char name[1024];
  snprintf(name, 1024, "JFNK Adams-Moulton (order %d)", order);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  am_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void* am_ode_integrator_context(ode_integrator_t* integrator)
{
  am_ode_t* integ = ode_integrator_context(integrator);
  return integ->context;
}

void am_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                      real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  am_ode_t* integ = ode_integrator_context(integrator);

  // Clear any existing error weight function.
  integ->compute_weights = NULL;

  // Set the tolerances.
  CVodeSStolerances(integ->cvode, relative_tol, absolute_tol);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  am_ode_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));
  return 0;
}

void am_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                 am_ode_integrator_error_weight_func compute_weights)
{
  am_ode_t* integ = ode_integrator_context(integrator);
  ASSERT(compute_weights != NULL);
  integ->compute_weights = compute_weights;
  CVodeWFtolerances(integ->cvode, compute_error_weights);
}

void am_ode_integrator_eval_rhs(ode_integrator_t* integrator, real_t t, real_t* X, real_t* rhs)
{
  am_ode_t* integ = ode_integrator_context(integrator);
  integ->rhs(integ->context, t, X, rhs);
}

void am_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                       am_ode_integrator_diagnostics_t* diagnostics)
{
  am_ode_t* integ = ode_integrator_context(integrator);
  diagnostics->status_message = integ->status_message; // borrowed!
  CVodeGetNumSteps(integ->cvode, &diagnostics->num_steps);
  CVodeGetLastOrder(integ->cvode, &diagnostics->order_of_last_step);
  CVodeGetLastStep(integ->cvode, &diagnostics->last_step_size);
  CVodeGetNumRhsEvals(integ->cvode, &diagnostics->num_rhs_evaluations);
  CVodeGetNumLinSolvSetups(integ->cvode, &diagnostics->num_linear_solve_setups);
  CVodeGetNumErrTestFails(integ->cvode, &diagnostics->num_error_test_failures);
  CVodeGetNumNonlinSolvIters(integ->cvode, &diagnostics->num_nonlinear_solve_iterations);
  CVodeGetNumNonlinSolvConvFails(integ->cvode, &diagnostics->num_nonlinear_solve_convergence_failures);
  CVSpilsGetNumLinIters(integ->cvode, &diagnostics->num_linear_solve_iterations);
  CVSpilsGetNumPrecEvals(integ->cvode, &diagnostics->num_preconditioner_evaluations);
  CVSpilsGetNumPrecSolves(integ->cvode, &diagnostics->num_preconditioner_solves);
  CVSpilsGetNumConvFails(integ->cvode, &diagnostics->num_linear_solve_convergence_failures);
}

void am_ode_integrator_diagnostics_fprintf(am_ode_integrator_diagnostics_t* diagnostics, 
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

