// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "integrators/nonlinear_integrator.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

// We use a serial version of SuperLU to do preconditioning.
#include "slu_ddefs.h"
#include "supermatrix.h"
#include "slu_util.h"

struct nonlinear_integrator_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  char* name;
  void* context;
  nonlinear_integrator_vtable vtable;
  nonlinear_integrator_type_t type;

  // KINSol data structures.
  void* kinsol;
  N_Vector x, x_scale, F_scale; // Stores solution vector and scaling vectors.

  // Preconditioning stuff.
  SuperMatrix precond_mat, precond_rhs, precond_L, precond_U;
  int *precond_rperm, *precond_cperm;
};

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(N_Vector x, N_Vector x_scale, 
                                 N_Vector F, N_Vector f_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  // FIXME: SuperLU stuff goes here!
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
  // FIXME: SuperLU stuff goes here!
  return 0;
}

nonlinear_integrator_t* nonlinear_integrator_new(const char* name, 
                                                 void* context,
                                                 MPI_Comm comm,
                                                 nonlinear_integrator_vtable vtable,
                                                 nonlinear_integrator_type_t type,
                                                 int max_krylov_dim)
{
  ASSERT(vtable.eval != NULL);
  ASSERT(vtable.graph != NULL);
  ASSERT(max_krylov_dim >= 3);

  nonlinear_integrator_t* integrator = malloc(sizeof(nonlinear_integrator_t));
  integrator->name = string_dup(name);
  integrator->context = context;
  integrator->comm = comm;
  integrator->vtable = vtable;
  integrator->type = type;

  // Get the adjacency graph that expresses the sparsity of the nonlinear system.
  adj_graph_t* graph = vtable.graph(context);
  ASSERT(graph != NULL);

  // The dimension N of the system is the number of local vertices in the 
  // adjacency graph.
  int N = adj_graph_num_vertices(graph);

  // Set up KINSol and accessories.
  integrator->x = N_VNew(comm, N);
  integrator->x_scale = N_VNew(comm, N);
  integrator->F_scale = N_VNew(comm, N);
  integrator->kinsol = KINCreate();
  KINSetUserData(integrator->kinsol, integrator->context);
  KINInit(integrator->kinsol, vtable.eval, integrator->x);

  // Set the constraints (if any) for the solution.
  if (integrator->vtable.set_constraints != NULL)
  {
    N_Vector constraints = N_VNew(comm, N);
    integrator->vtable.set_constraints(integrator->context, NV_DATA(constraints));
    KINSetConstraints(integrator->kinsol, constraints);
    N_VDestroy(constraints);
  }

  // Select the particular type of Krylov method for the underlying linear solves.
  if (type == GMRES)
    KINSpgmr(integrator->kinsol, max_krylov_dim); 
  else if (type == BICGSTAB)
    KINSpbcg(integrator->kinsol, max_krylov_dim);
  else
    KINSptfqmr(integrator->kinsol, max_krylov_dim);

  // Set up preconditioner machinery.
  KINSpilsSetPreconditioner(integrator->kinsol, set_up_preconditioner,
                            solve_preconditioner_system);

  return integrator;
}

// This function frees all of the machinery having to do with the preconditioner.
static void free_preconditioner(nonlinear_integrator_t* integrator)
{
  SUPERLU_FREE(integrator->precond_cperm);
  SUPERLU_FREE(integrator->precond_rperm);
//  Destroy_CompCol_Matrix(&mf_integrator->precond_U);
//  Destroy_SuperNode_Matrix(&mf_integrator->precond_L);
  Destroy_SuperMatrix_Store(&integrator->precond_rhs);
  Destroy_CompRow_Matrix(&integrator->precond_mat);
}

void nonlinear_integrator_free(nonlinear_integrator_t* integrator)
{
//  free_preconditioner(integrator);
  N_VDestroy(integrator->x);
  N_VDestroy(integrator->x_scale);
  N_VDestroy(integrator->F_scale);
  KINFree(&integrator->kinsol);
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

void nonlinear_integrator_set_tolerances(nonlinear_integrator_t* integrator, double norm_tolerance, double step_tolerance)
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

bool nonlinear_integrator_solve(nonlinear_integrator_t* integrator,
                                double t,
                                double* X,
                                int* num_iterations)
{
  ASSERT(X != NULL);

  // Get the adjacency graph that expresses the sparsity of the nonlinear system.
  adj_graph_t* graph = integrator->vtable.graph(integrator->context);
  ASSERT(graph != NULL);

  // The dimension N of the system is the number of local vertices in the 
  // adjacency graph.
  int N = adj_graph_num_vertices(graph);
  ASSERT(NV_LOCLENGTH(integrator->x) == N); // No adaptivity allowed yet!

  // Set the current time in the state.
  if (integrator->vtable.set_time != NULL)
    integrator->vtable.set_time(integrator->context, t);

  // Set the x_scale and F_scale vectors. If we don't have methods for doing 
  // this, the scaling vectors are set to 1.
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
  memcpy(NV_DATA(integrator->x), X, sizeof(double) * N);

  // Suspend the currently active floating point exceptions for now.
  polymec_suspend_fpe_exceptions();

  // Solve.
  int status = KINSol(integrator->kinsol, integrator->x, KIN_LINESEARCH, 
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
    memcpy(X, NV_DATA(integrator->x), sizeof(double) * N);
    return true;
  }

  // Failed!
  return false;
}
                                  
