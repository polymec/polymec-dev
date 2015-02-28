// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/sundials_helpers.h"
#include "integrators/krylov_solver.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

struct krylov_solver_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  void* context;
  krylov_t solver_type;
  int max_krylov_dim, max_restarts;

  int num_local_values, num_remote_values;

  int (*Ay)(void* context, real_t* y, real_t* Ay);
  void (*dtor)(void* context);

  // KINSol data structures.
  void* kinsol;
  N_Vector x, x_scale, F_scale; // Stores solution vector, scaled vectors.
  real_t* x_with_ghosts;
  real_t* b;
  char* status_message; // status of most recent integration.

  // Preconditioning stuff.
//  newton_pc_t* precond;
};

// This function wraps around the user-supplied evaluation function.
static int evaluate_R(N_Vector x, N_Vector R, void* context)
{
  krylov_solver_t* solver = context;
  real_t* xx = NV_DATA(x);
  real_t* RR = NV_DATA(R);

  // Evaluate the residual using a solution vector with ghosts.
  memcpy(solver->x_with_ghosts, xx, sizeof(real_t) * solver->num_local_values);
  int status = solver->Ay(solver->context, solver->x_with_ghosts, RR);
  if (status == 0)
  {
    // Subtract b.
    for (int i = 0; i < solver->num_local_values; ++i)
      RR[i] -= solver->b[i];
  }

  return status;
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(N_Vector x, N_Vector x_scale, 
                                 N_Vector F, N_Vector f_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  krylov_solver_t* solver = context;
  log_debug("krylov_solver: setting up preconditioner...");
//  newton_pc_setup(solver->precond, 0.0, 1.0, 0.0, t, NV_DATA(x), NULL);

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
  krylov_solver_t* solver = context;

  // FIXME: Apply scaling if needed.

  log_debug("krylov_solver: solving preconditioner...");
//  if (newton_pc_solve(solver->precond, NV_DATA(r)))
    return 0;
#if 0
  else 
  {
    // Recoverable error.
    log_debug("krylov_solver: preconditioner solve failed.");
    return 1; 
  }
#endif
}

// Generic constructor.
krylov_solver_t* krylov_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   int (*matrix_vector_product)(void* context, real_t* y, real_t* Ay),
                                   void (*dtor)(void* context),
                                   krylov_t solver_type,
                                   int max_krylov_dim,
                                   int max_restarts)
{
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(matrix_vector_product != NULL);
  ASSERT(max_krylov_dim >= 3);
  ASSERT((solver_type != KRYLOV_GMRES) || (max_restarts >= 0));

  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->context = context;
  solver->comm = comm;
  solver->Ay = matrix_vector_product;
  solver->dtor = dtor;
//  solver->precond = precond;
  solver->solver_type = solver_type;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->max_krylov_dim = max_krylov_dim;
  solver->max_restarts = max_restarts;

  // Set up KINSol and accessories.
  solver->kinsol = KINCreate();
  KINSetUserData(solver->kinsol, solver);
  solver->x = N_VNew(solver->comm, num_local_values);
  solver->x_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->b = polymec_malloc(sizeof(real_t) * solver->num_local_values);
  solver->status_message = NULL;

  KINInit(solver->kinsol, evaluate_R, solver->x);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver->solver_type == KRYLOV_GMRES)
  {
    KINSpgmr(solver->kinsol, solver->max_krylov_dim); 
    KINSpilsSetMaxRestarts(solver->kinsol, solver->max_restarts);
  }
  else if (solver->solver_type == KRYLOV_BICGSTAB)
    KINSpbcg(solver->kinsol, solver->max_krylov_dim);
  else
    KINSptfqmr(solver->kinsol, solver->max_krylov_dim);

  // Enable debugging diagnostics if logging permits.
  FILE* info_stream = log_stream(LOG_DEBUG);
  if (info_stream != NULL)
  {
    KINSetPrintLevel(solver->kinsol, 3);
    KINSetInfoFile(solver->kinsol, info_stream);
  }
  else
  {
    KINSetPrintLevel(solver->kinsol, 0);
    KINSetInfoFile(solver->kinsol, NULL);
  }

  // Set up the preconditioner.
//  if (solver->precond != NULL)
//  {
//    KINSpilsSetPreconditioner(solver->kinsol, set_up_preconditioner,
//                              solve_preconditioner_system);
//  }

  return solver;
}

void krylov_solver_free(krylov_solver_t* solver)
{
  // Kill the preconditioner stuff.
//  if (solver->precond != NULL)
//    newton_pc_free(solver->precond);

  // Kill the KINSol stuff.
  N_VDestroy(solver->x);
  KINFree(&solver->kinsol);

  // Kill the rest.
  if ((solver->dtor != NULL) && (solver->context != NULL))
    solver->dtor(solver->context);
  // Kill the rest.
  if (solver->status_message != NULL)
    polymec_free(solver->status_message);
  polymec_free(solver->x_with_ghosts);
  polymec_free(solver->b);
  polymec_free(solver);
}

void* krylov_solver_context(krylov_solver_t* solver)
{
  return solver->context;
}

int krylov_solver_num_equations(krylov_solver_t* solver)
{
  return solver->num_local_values;
}

void krylov_solver_set_tolerances(krylov_solver_t* solver, real_t norm_tolerance, real_t step_tolerance)
{
  ASSERT(norm_tolerance > 0.0);
  ASSERT(step_tolerance > 0.0);
  KINSetFuncNormTol(solver->kinsol, norm_tolerance);
  KINSetScaledStepTol(solver->kinsol, step_tolerance);
}

void krylov_solver_set_max_iterations(krylov_solver_t* solver, int max_iterations)
{
  ASSERT(max_iterations > 0);
  KINSetNumMaxIters(solver->kinsol, max_iterations);
}

bool krylov_solver_solve(krylov_solver_t* solver,
                         real_t* b,
                         int* num_iterations)
{
  ASSERT(b != NULL);

  // Copy b into place.
  int N = solver->num_local_values;
  memcpy(solver->b, b, sizeof(real_t) * N);

  // Set the x_scale and F_scale vectors. If we don't have methods for doing 
  // this, the scaling vectors are set to 1.
  for (int i = 0; i < N; ++i)
    NV_Ith(solver->x_scale, i) = 1.0;

  for (int i = 0; i < N; ++i)
    NV_Ith(solver->F_scale, i) = 1.0;

  // Zero the internal solution vector.
  memset(NV_DATA(solver->x), 0, sizeof(real_t) * N);

  // Suspend the currently active floating point exceptions for now.
//  polymec_suspend_fpe_exceptions();

  // Solve.
  log_debug("krylov_solver: solving...");
  int status = KINSol(solver->kinsol, solver->x, KIN_NONE, solver->x_scale, solver->F_scale);

  // Clear the present status.
  if (solver->status_message != NULL)
  {
    polymec_free(solver->status_message);
    solver->status_message = NULL;
  }

  // Reinstate the floating point exceptions.
//  polymec_restore_fpe_exceptions();

  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(solver->kinsol, &num_iters);
    *num_iterations = (int)num_iters;
    log_debug("krylov_solver: solved after %d iterations.", *num_iterations);

    // Copy the data back into b.
    memcpy(b, NV_DATA(solver->x), sizeof(real_t) * N);
    return true;
  }
  else
  {
    if (status == KIN_STEP_LT_STPTOL)
      solver->status_message = string_dup("Nonlinear solve stalled because scaled Newton step is too small.");
    else if (status == KIN_LINESEARCH_NONCONV)
      solver->status_message = string_dup("Line search could not sufficiently decrease the error of the iterate.");
    else if (status == KIN_MAXITER_REACHED)
      solver->status_message = string_dup("Maximum number of nonlinear iterations was reached.");
    else if (status == KIN_MXNEWT_5X_EXCEEDED)
      solver->status_message = string_dup("Maximum Newton step size was exceeded 5 times.");
    else if (status == KIN_LINESEARCH_BCFAIL)
      solver->status_message = string_dup("Line search could not satisfy beta condition.");
    else if (status == KIN_LINSOLV_NO_RECOVERY)
      solver->status_message = string_dup("Preconditioner solve encountered a recoverable error after update.");
    else if (status == KIN_LINIT_FAIL)
      solver->status_message = string_dup("Linear solve setup failed.");
    else if (status == KIN_LSETUP_FAIL)
      solver->status_message = string_dup("Preconditioner setup failed unrecoverably.");
    else if (status == KIN_LSOLVE_FAIL)
      solver->status_message = string_dup("Linear solve failed (or preconditioner solve failed unrecoverably).");
    else if (status == KIN_SYSFUNC_FAIL)
      solver->status_message = string_dup("Nonlinear function evaluation failed unrecoverably.");
    else if (status == KIN_FIRST_SYSFUNC_ERR)
      solver->status_message = string_dup("First nonlinear function evaluation failed recoverably.");
    else if (status == KIN_REPTD_SYSFUNC_ERR)
      solver->status_message = string_dup("Nonlinear function evaluation repeatedly failed (no recovery possible).");
  }

  // Failed!
  return false;
}
                                  
void krylov_solver_get_diagnostics(krylov_solver_t* solver, 
                                   krylov_solver_diagnostics_t* diagnostics)
{
  diagnostics->status_message = solver->status_message; // borrowed!
  KINGetNumFuncEvals(solver->kinsol, &diagnostics->num_function_evaluations);
  KINGetFuncNorm(solver->kinsol, &diagnostics->function_norm);
  KINSpilsGetNumLinIters(solver->kinsol, &diagnostics->num_solve_iterations);
  KINSpilsGetNumConvFails(solver->kinsol, &diagnostics->num_solve_convergence_failures);
  KINSpilsGetNumPrecEvals(solver->kinsol, &diagnostics->num_preconditioner_evaluations);
  KINSpilsGetNumPrecSolves(solver->kinsol, &diagnostics->num_preconditioner_solves);
  KINSpilsGetNumJtimesEvals(solver->kinsol, &diagnostics->num_matrix_vector_product_evaluations);
}

void krylov_solver_diagnostics_fprintf(krylov_solver_diagnostics_t* diagnostics, 
                                       FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "Nonlinear solver diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num function evaluations: %d\n", (int)diagnostics->num_function_evaluations);
  fprintf(stream, "  Function norm: %g\n", (double)diagnostics->function_norm);
  fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_solve_iterations);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_solve_convergence_failures);
  fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
  fprintf(stream, "  Num matrix-vector product evaluations: %d\n", (int)diagnostics->num_matrix_vector_product_evaluations);
}

