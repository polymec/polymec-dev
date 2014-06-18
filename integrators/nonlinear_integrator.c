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

#include <float.h>
#include "core/sundials_helpers.h"
#include "integrators/nonlinear_integrator.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

// Types of linear solvers.
typedef enum 
{
  GMRES,
  BICGSTAB,
  TFQMR
} solver_type_t;

struct nonlinear_integrator_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  char* name;
  void* context;
  nonlinear_integrator_vtable vtable;
  solver_type_t solver_type;
  int max_krylov_dim, max_restarts;

  int N; // Number of degrees of freedom.

  // KINSol data structures.
  void* kinsol;
  int strategy; // Global strategy.
  N_Vector x, x_scale, F_scale; // Stores solution vector and scaling vectors.
  char* status_message; // status of most recent integration.

  // Preconditioning stuff.
  preconditioner_t* precond;

  // Null space information.
  bool homogeneous_functions_in_null_space;
  real_t** null_space_vectors;
  int null_dim;

  // Current simulation time.
  real_t current_time;
};

static void project_out_of_null_space(nonlinear_integrator_t* integrator,
                                      real_t* R)
{
  // If homogeneous functions are in the null space, subtract the spatial 
  // mean from R.
  if (integrator->homogeneous_functions_in_null_space)
  {
    real_t mean = 0.0;
    for (int i = 0; i < integrator->N; ++i)
      mean += R[i];
    mean *= 1.0/integrator->N;
    for (int i = 0; i < integrator->N; ++i)
      R[i] -= mean;
  }

  // FIXME: Do the other stuff here.
}

// This function wraps around the user-supplied evaluation function.
static int evaluate_F(N_Vector x, N_Vector F, void* context)
{
  nonlinear_integrator_t* integrator = context;
  real_t* xx = NV_DATA(x);
  real_t* FF = NV_DATA(F);

  // Evaluate the residual.
  int retval = integrator->vtable.eval(integrator->context, integrator->current_time, xx, FF);
  project_out_of_null_space(integrator, FF);
  return retval;
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(N_Vector x, N_Vector x_scale, 
                                 N_Vector F, N_Vector f_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  nonlinear_integrator_t* integrator = context;
  real_t t = integrator->current_time;
  preconditioner_setup(integrator->precond, 0.0, 1.0, 0.0, t, NV_DATA(x), NULL);
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(N_Vector x, N_Vector x_scale,
                                       N_Vector F, N_Vector F_scale,
                                       N_Vector r, void* context,
                                       N_Vector work)
{
  nonlinear_integrator_t* integrator = context;

  // FIXME: Apply scaling if needed.

  // Project r out of the null space.
  project_out_of_null_space(integrator, NV_DATA(r));

  if (preconditioner_solve(integrator->precond, NV_DATA(r)))
    return 0;
  else 
  {
    // Recoverable error.
    log_debug("nonlinear_integrator: preconditioner solve failed.");
    return 1; 
  }
}

// Generic constructor.
static nonlinear_integrator_t* nonlinear_integrator_new(const char* name, 
                                                        void* context,
                                                        MPI_Comm comm,
                                                        int N,
                                                        nonlinear_integrator_vtable vtable,
                                                        nonlinear_integrator_strategy_t global_strategy,
                                                        solver_type_t solver_type,
                                                        int max_krylov_dim, 
                                                        int max_restarts)
{
  ASSERT(N > 0);
  ASSERT(vtable.eval != NULL);
  ASSERT(max_krylov_dim >= 3);
  ASSERT(max_restarts >= 0);

  nonlinear_integrator_t* integrator = polymec_malloc(sizeof(nonlinear_integrator_t));
  integrator->name = string_dup(name);
  integrator->context = context;
  integrator->comm = comm;
  integrator->vtable = vtable;
  integrator->solver_type = solver_type;
  integrator->strategy = (global_strategy == LINE_SEARCH) ? KIN_LINESEARCH : KIN_NONE;
  integrator->N = N;
  integrator->max_krylov_dim = max_krylov_dim;
  integrator->max_restarts = max_restarts;

  // Set up KINSol and accessories.
  integrator->kinsol = KINCreate();
  KINSetUserData(integrator->kinsol, integrator);
  integrator->x = N_VNew(integrator->comm, N);
  integrator->x_scale = N_VNew(integrator->comm, N);
  integrator->F_scale = N_VNew(integrator->comm, N);
  integrator->status_message = NULL;

  KINInit(integrator->kinsol, evaluate_F, integrator->x);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (integrator->solver_type == GMRES)
  {
    KINSpgmr(integrator->kinsol, integrator->max_krylov_dim); 
    KINSpilsSetMaxRestarts(integrator->kinsol, integrator->max_restarts);
  }
  else if (integrator->solver_type == BICGSTAB)
    KINSpbcg(integrator->kinsol, integrator->max_krylov_dim);
  else
    KINSptfqmr(integrator->kinsol, integrator->max_krylov_dim);

  // Set up the preconditioner.
  KINSpilsSetPreconditioner(integrator->kinsol, set_up_preconditioner,
                            solve_preconditioner_system);

  // Enable debugging diagnostics if logging permits.
  FILE* info_stream = log_stream(LOG_DEBUG);
  if (info_stream != NULL)
  {
    KINSetPrintLevel(integrator->kinsol, 3);
    KINSetInfoFile(integrator->kinsol, info_stream);
  }
  else
  {
    KINSetPrintLevel(integrator->kinsol, 0);
    KINSetInfoFile(integrator->kinsol, NULL);
  }

  // Set the constraints (if any) for the solution.
  if (integrator->vtable.set_constraints != NULL)
  {
    N_Vector constraints = N_VNew(integrator->comm, N);
    integrator->vtable.set_constraints(integrator->context, NV_DATA(constraints));
    KINSetConstraints(integrator->kinsol, constraints);
    N_VDestroy(constraints);
  }

  integrator->precond = NULL;
  integrator->current_time = 0.0;

  // Set up the null space.
  integrator->homogeneous_functions_in_null_space = false;
  integrator->null_space_vectors = NULL;
  integrator->null_dim = 0;

  return integrator;
}

nonlinear_integrator_t* gmres_nonlinear_integrator_new(const char* name,
                                                       void* context,
                                                       MPI_Comm comm,
                                                       int N,
                                                       nonlinear_integrator_vtable vtable,
                                                       nonlinear_integrator_strategy_t global_strategy,
                                                       int max_krylov_dim,
                                                       int max_restarts)
{
  return nonlinear_integrator_new(name, context, comm, N, vtable, global_strategy,
                                  GMRES, max_krylov_dim, max_restarts);
}

nonlinear_integrator_t* bicgstab_nonlinear_integrator_new(const char* name,
                                                          void* context,
                                                          MPI_Comm comm,
                                                          int N,
                                                          nonlinear_integrator_vtable vtable,
                                                          nonlinear_integrator_strategy_t global_strategy,
                                                          int max_krylov_dim)
{
  return nonlinear_integrator_new(name, context, comm, N, vtable, global_strategy,
                                  BICGSTAB, max_krylov_dim, 0);
}

nonlinear_integrator_t* tfqmr_nonlinear_integrator_new(const char* name,
                                                       void* context,
                                                       MPI_Comm comm,
                                                       int N,
                                                       nonlinear_integrator_vtable vtable,
                                                       nonlinear_integrator_strategy_t global_strategy,
                                                       int max_krylov_dim)
{
  return nonlinear_integrator_new(name, context, comm, N, vtable, global_strategy,
                                  TFQMR, max_krylov_dim, 0);
}

void nonlinear_integrator_free(nonlinear_integrator_t* integrator)
{
  // Kill the null space.
  nonlinear_integrator_set_null_space(integrator, false, NULL, 0);

  // Kill the preconditioner stuff.
  if (integrator->precond != NULL)
    preconditioner_free(integrator->precond);

  // Kill the KINSol stuff.
  N_VDestroy(integrator->x);
  N_VDestroy(integrator->x_scale);
  N_VDestroy(integrator->F_scale);
  KINFree(&integrator->kinsol);

  // Kill the rest.
  if ((integrator->vtable.dtor != NULL) && (integrator->context != NULL))
    integrator->vtable.dtor(integrator->context);
  // Kill the rest.
  if (integrator->status_message != NULL)
    polymec_free(integrator->status_message);
  polymec_free(integrator->name);
  polymec_free(integrator);
}

char* nonlinear_integrator_name(nonlinear_integrator_t* integrator)
{
  return integrator->name;
}

void* nonlinear_integrator_context(nonlinear_integrator_t* integrator)
{
  return integrator->context;
}

int nonlinear_integrator_num_equations(nonlinear_integrator_t* integrator)
{
  return integrator->N;
}

void nonlinear_integrator_set_tolerances(nonlinear_integrator_t* integrator, real_t norm_tolerance, real_t step_tolerance)
{
  ASSERT(norm_tolerance > 0.0);
  ASSERT(step_tolerance > 0.0);
  KINSetFuncNormTol(integrator->kinsol, norm_tolerance);
  KINSetScaledStepTol(integrator->kinsol, step_tolerance);
}

void newton_solver_set_max_iterations(nonlinear_integrator_t* integrator, int max_iterations)
{
  ASSERT(max_iterations > 0);
  KINSetNumMaxIters(integrator->kinsol, max_iterations);
}

void nonlinear_integrator_set_preconditioner(nonlinear_integrator_t* integrator,
                                             preconditioner_t* precond)
{
  integrator->precond = precond;
}

preconditioner_t* nonlinear_integrator_preconditioner(nonlinear_integrator_t* integrator)
{
  return integrator->precond;
}

void nonlinear_integrator_set_null_space(nonlinear_integrator_t* integrator,
                                         bool homogeneous_functions,
                                         real_t** null_space_vectors,
                                         int null_dim)
{
  ASSERT(((null_dim == 0) && (null_space_vectors == NULL)) ||
         ((null_dim > 0) && (null_space_vectors != NULL)));

  integrator->homogeneous_functions_in_null_space = homogeneous_functions;
  if (integrator->null_space_vectors != NULL)
  {
    for (int i = 0; i < integrator->null_dim; ++i)
      polymec_free(integrator->null_space_vectors[i]);
    polymec_free(integrator->null_space_vectors);
    integrator->null_dim = 0;
    integrator->null_space_vectors = NULL;
  }
  if (null_space_vectors != NULL)
  {
    integrator->null_dim = null_dim;
    integrator->null_space_vectors = polymec_malloc(sizeof(real_t*));
    for (int i = 0; i < null_dim; ++i)
    {
      integrator->null_space_vectors[i] = polymec_malloc(sizeof(real_t) * integrator->N);
      memcpy(integrator->null_space_vectors[i], null_space_vectors[i], integrator->N * sizeof(real_t));
    }
  }
}

void nonlinear_integrator_eval_residual(nonlinear_integrator_t* integrator, real_t t, real_t* X, real_t* F)
{
  integrator->vtable.eval(integrator->context, t, X, F);
  project_out_of_null_space(integrator, F);
}

bool nonlinear_integrator_solve(nonlinear_integrator_t* integrator,
                                real_t t,
                                real_t* X,
                                int* num_iterations)
{
  ASSERT(X != NULL);

  // Set the current time in the state.
  integrator->current_time = t;

  // Set the x_scale and F_scale vectors. If we don't have methods for doing 
  // this, the scaling vectors are set to 1.
  int N = integrator->N;
  if (integrator->vtable.set_x_scale != NULL)
    integrator->vtable.set_x_scale(integrator->context, NV_DATA(integrator->x_scale));
  else
  {
    for (int i = 0; i < N; ++i)
      NV_Ith(integrator->x_scale, i) = 1.0;
  }

  if (integrator->vtable.set_F_scale != NULL)
    integrator->vtable.set_F_scale(integrator->context, NV_DATA(integrator->F_scale));
  else
  {
    for (int i = 0; i < N; ++i)
      NV_Ith(integrator->F_scale, i) = 1.0;
  }

  if (integrator->vtable.initial_guess != NULL)
  {
    // Form the initial guess magically.
    log_debug("nonlinear_integrator: forming initial guess...");
    integrator->vtable.initial_guess(integrator->context, t, NV_DATA(integrator->x));
  }
  else
  {
    // Copy the values in X to the internal solution vector.
    memcpy(NV_DATA(integrator->x), X, sizeof(real_t) * N);
  }

  // Suspend the currently active floating point exceptions for now.
//  polymec_suspend_fpe_exceptions();

  // Solve.
  log_debug("nonlinear_integrator: solving...");
  int status = KINSol(integrator->kinsol, integrator->x, integrator->strategy, 
                      integrator->x_scale, integrator->F_scale);

  // Clear the present status.
  if (integrator->status_message != NULL)
  {
    polymec_free(integrator->status_message);
    integrator->status_message = NULL;
  }

  // Reinstate the floating point exceptions.
//  polymec_restore_fpe_exceptions();

  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(integrator->kinsol, &num_iters);
    *num_iterations = (int)num_iters;
    log_debug("nonlinear_integrator: solved after %d iterations.", *num_iterations);

    // Copy the data back into X.
    memcpy(X, NV_DATA(integrator->x), sizeof(real_t) * N);
    return true;
  }
  else
  {
    if (status == KIN_STEP_LT_STPTOL)
      integrator->status_message = string_dup("Nonlinear solve stalled because scaled Newton step is too small.");
    else if (status == KIN_LINESEARCH_NONCONV)
      integrator->status_message = string_dup("Line search could not sufficiently decrease the error of the iterate.");
    else if (status == KIN_MAXITER_REACHED)
      integrator->status_message = string_dup("Maximum number of nonlinear iterations was reached.");
    else if (status == KIN_MXNEWT_5X_EXCEEDED)
      integrator->status_message = string_dup("Maximum Newton step size was exceeded 5 times.");
    else if (status == KIN_LINESEARCH_BCFAIL)
      integrator->status_message = string_dup("Line search could not satisfy beta condition.");
    else if (status == KIN_LINSOLV_NO_RECOVERY)
      integrator->status_message = string_dup("Preconditioner solve encountered a recoverable error after update.");
    else if (status == KIN_LINIT_FAIL)
      integrator->status_message = string_dup("Linear solve setup failed.");
    else if (status == KIN_LSETUP_FAIL)
      integrator->status_message = string_dup("Preconditioner setup failed unrecoverably.");
    else if (status == KIN_LSOLVE_FAIL)
      integrator->status_message = string_dup("Linear solve failed (or preconditioner solve failed unrecoverably).");
    else if (status == KIN_SYSFUNC_FAIL)
      integrator->status_message = string_dup("Nonlinear function evaluation failed unrecoverably.");
    else if (status == KIN_FIRST_SYSFUNC_ERR)
      integrator->status_message = string_dup("First nonlinear function evaluation failed recoverably.");
    else if (status == KIN_REPTD_SYSFUNC_ERR)
      integrator->status_message = string_dup("Nonlinear function evaluation repeatedly failed (no recovery possible).");
  }

  // Failed!
  return false;
}
                                  
void nonlinear_integrator_get_diagnostics(nonlinear_integrator_t* integrator, 
                                          nonlinear_integrator_diagnostics_t* diagnostics)
{
  diagnostics->status_message = integrator->status_message; // borrowed!
  KINGetNumFuncEvals(integrator->kinsol, &diagnostics->num_function_evaluations);
  KINGetNumBetaCondFails(integrator->kinsol, &diagnostics->num_beta_condition_failures);
  KINGetNumNonlinSolvIters(integrator->kinsol, &diagnostics->num_nonlinear_iterations);
  KINGetNumBacktrackOps(integrator->kinsol, &diagnostics->num_backtrack_operations);
  KINGetFuncNorm(integrator->kinsol, &diagnostics->scaled_function_norm);
  KINGetStepLength(integrator->kinsol, &diagnostics->scaled_newton_step_length);
  KINSpilsGetNumLinIters(integrator->kinsol, &diagnostics->num_linear_solve_iterations);
  KINSpilsGetNumConvFails(integrator->kinsol, &diagnostics->num_linear_solve_convergence_failures);
  KINSpilsGetNumPrecEvals(integrator->kinsol, &diagnostics->num_preconditioner_evaluations);
  KINSpilsGetNumPrecSolves(integrator->kinsol, &diagnostics->num_preconditioner_solves);
  KINSpilsGetNumJtimesEvals(integrator->kinsol, &diagnostics->num_jacobian_vector_product_evaluations);
  KINSpilsGetNumFuncEvals(integrator->kinsol, &diagnostics->num_difference_quotient_function_evaluations);
}

void nonlinear_integrator_diagnostics_fprintf(nonlinear_integrator_diagnostics_t* diagnostics, 
                                              FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "Nonlinear integrator diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num function evaluations: %d\n", (int)diagnostics->num_function_evaluations);
  fprintf(stream, "  Num beta condition failures: %d\n", (int)diagnostics->num_beta_condition_failures);
  fprintf(stream, "  Num backtrack operations: %d\n", (int)diagnostics->num_backtrack_operations);
  fprintf(stream, "  Num nonlinear iterations: %d\n", (int)diagnostics->num_nonlinear_iterations);
  fprintf(stream, "  Scaled function norm: %g\n", (double)diagnostics->scaled_function_norm);
  fprintf(stream, "  Scaled Newton step length: %g\n", (double)diagnostics->scaled_newton_step_length);
  fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
  fprintf(stream, "  Num Jacobian-vector product evaluations: %d\n", (int)diagnostics->num_jacobian_vector_product_evaluations);
  fprintf(stream, "  Num difference quotient function evaluations: %d\n", (int)diagnostics->num_difference_quotient_function_evaluations);
}

