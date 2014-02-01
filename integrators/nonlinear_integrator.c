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

  // Preconditioning stuff.
  preconditioner_t* precond;
  preconditioner_matrix_t* precond_mat;

  // Current simulation time.
  real_t current_time;
};

// This function wraps around the user-supplied evaluation function.
static int evaluate_F(N_Vector x, N_Vector F, void* context)
{
  nonlinear_integrator_t* integrator = context;
  real_t* xx = NV_DATA(x);
  real_t* FF = NV_DATA(F);

  // Do parallel communication.
  if (integrator->vtable.communicate != NULL)
    integrator->vtable.communicate(integrator->context, integrator->current_time, xx);

  // Evaluate the residual.
  return integrator->vtable.eval(integrator->context, integrator->current_time, xx, FF);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(N_Vector x, N_Vector x_scale, 
                                 N_Vector F, N_Vector f_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  nonlinear_integrator_t* integrator = context;
  real_t t = integrator->current_time;
  preconditioner_compute_jacobian(integrator->precond, t, NV_DATA(x), integrator->precond_mat);
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

  preconditioner_solve(integrator->precond, integrator->precond_mat, NV_DATA(r));

  return 0;
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

  nonlinear_integrator_t* integrator = malloc(sizeof(nonlinear_integrator_t));
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

  KINSpilsSetPreconditioner(integrator->kinsol, set_up_preconditioner,
                            solve_preconditioner_system);

  // Set the constraints (if any) for the solution.
  if (integrator->vtable.set_constraints != NULL)
  {
    N_Vector constraints = N_VNew(integrator->comm, N);
    integrator->vtable.set_constraints(integrator->context, NV_DATA(constraints));
    KINSetConstraints(integrator->kinsol, constraints);
    N_VDestroy(constraints);
  }

  integrator->precond = NULL;
  integrator->precond_mat = NULL;
  integrator->current_time = 0.0;

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
  // Kill the preconditioner stuff.
  if (integrator->precond != NULL)
    preconditioner_free(integrator->precond);
  if (integrator->precond_mat != NULL)
    preconditioner_matrix_free(integrator->precond_mat);

  // Kill the KINSol stuff.
  N_VDestroy(integrator->x);
  N_VDestroy(integrator->x_scale);
  N_VDestroy(integrator->F_scale);
  KINFree(&integrator->kinsol);

  // Kill the rest.
  if ((integrator->vtable.dtor != NULL) && (integrator->context != NULL))
    integrator->vtable.dtor(integrator->context);
  free(integrator->name);
  free(integrator);
}

char* nonlinear_integrator_name(nonlinear_integrator_t* integrator)
{
  return integrator->name;
}

void* nonlinear_integrator_context(nonlinear_integrator_t* integrator)
{
  return integrator->context;
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
  if (integrator->precond_mat != NULL)
    preconditioner_matrix_free(integrator->precond_mat);
  integrator->precond_mat = preconditioner_matrix(precond);
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

  // Copy the values in X to the internal solution vector.
  memcpy(NV_DATA(integrator->x), X, sizeof(real_t) * N);

  // Suspend the currently active floating point exceptions for now.
  polymec_suspend_fpe_exceptions();

  // Solve.
  int status = KINSol(integrator->kinsol, integrator->x, integrator->strategy, 
                      integrator->x_scale, integrator->F_scale);

  // Reinstate the floating point exceptions.
  polymec_restore_fpe_exceptions();

  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(integrator->kinsol, &num_iters);
    *num_iterations = (int)num_iters;

    // Copy the data back into X.
    memcpy(X, NV_DATA(integrator->x), sizeof(real_t) * N);
    return true;
  }

  // Failed!
  return false;
}
                                  
