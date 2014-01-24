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
#include "integrators/time_integrator.h"
#include "integrators/supermatrix_factory.h"

#include "cvode/cvode.h"
#include "cvode/cvode_spils.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"
#include "cvode/cvode_sptfqmr.h"

struct time_integrator_t 
{
  char* name;
  void* context;
  time_integrator_vtable vtable;
  int order;
  MPI_Comm comm;
  time_integrator_solver_type_t solver_type;

  // Adjacency graph -- keeps track of topological changes.
  adj_graph_t* graph;

  // CVODE data structures.
  void* cvode;
  N_Vector x; 
  real_t current_time;

  // Preconditioning stuff.
  supermatrix_factory_t* precond_factory;
  SuperMatrix *precond_mat, *precond_rhs, precond_L, precond_U;
  int *precond_rperm, *precond_cperm;
  superlu_options_t precond_options;
  SuperLUStat_t precond_stat;
};

// This function wraps around the user-supplied right hand side.
static int evaluate_rhs(real_t t, N_Vector x, N_Vector x_dot, void* context)
{
  time_integrator_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);
  return integ->vtable.rhs(integ->context, t, xx, xxd);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  time_integrator_t* integ = context;
  if (!jacobian_is_current)
  {
    supermatrix_factory_update_jacobian(integ->precond_factory, 
                                        NV_DATA(x), t, integ->precond_mat);
    *jacobian_was_updated = 1;
    // FIXME: Incorporate gamma
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
  time_integrator_t* integ = context;
  int N = NV_LOCLENGTH(x); // Dimension of the matrix.
  
  // Copy the values from the vector r to the preconditioner right-hand side.
  {
    real_t *rhs = (real_t*) ((DNformat*) integ->precond_rhs->Store)->nzval; 
    memcpy(rhs, NV_DATA(x), sizeof(real_t) * N);
  }

  // Solve the preconditioner system.
  int info;
  dgssv(&integ->precond_options, integ->precond_mat, integ->precond_cperm,
        integ->precond_rperm, &integ->precond_L, 
        &integ->precond_U, integ->precond_rhs,
        &integ->precond_stat, &info);

  // Tell SuperLU to use the same nonzero pattern for the next factorization.
  integ->precond_options.Fact = SamePattern;

  // Copy the values from the solution to the vector r.
  {
    real_t *sol = (real_t*) ((DNformat*) integ->precond_rhs->Store)->nzval; 
    memcpy(NV_DATA(x), sol, sizeof(real_t) * N);
  }

  return 0;
}

time_integrator_t* time_integrator_new(const char* name, 
                                       void* context,
                                       MPI_Comm comm,
                                       time_integrator_vtable vtable,
                                       int order,
                                       time_integrator_solver_type_t solver_type,
                                       int max_krylov_dim)
{
  ASSERT(order > 0);
  ASSERT(vtable.rhs != NULL);
  ASSERT(vtable.graph != NULL);
  ASSERT(max_krylov_dim >= 3);

  time_integrator_t* integ = malloc(sizeof(time_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->comm = comm;
  integ->vtable = vtable;
  integ->order = order;
  integ->solver_type = solver_type;
  integ->graph = NULL;
  integ->current_time = 0.0;

  // Get the adjacency graph that expresses the sparsity of the nonlinear system.
  adj_graph_t* graph = vtable.graph(context);
  ASSERT(graph != NULL);

  // The dimension N of the system is the number of local vertices in the 
  // adjacency graph.
  int N = adj_graph_num_vertices(graph);

  // Set up KINSol and accessories.
  integ->x = N_VNew(comm, N);
  integ->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, evaluate_rhs, integ->current_time, integ->x);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver_type == GMRES)
    CVSpgmr(integ->cvode, PREC_LEFT, max_krylov_dim); 
  else if (solver_type == BICGSTAB)
    CVSpbcg(integ->cvode, PREC_LEFT, max_krylov_dim);
  else
    CVSptfqmr(integ->cvode, PREC_LEFT, max_krylov_dim);

  // Set up preconditioner machinery.
  CVSpilsSetPreconditioner(integ->cvode, set_up_preconditioner,
                           solve_preconditioner_system);
  set_default_options(&integ->precond_options);
  StatInit(&integ->precond_stat);
  integ->precond_factory = NULL;
  integ->precond_mat = NULL;
  integ->precond_rhs = NULL;

  return integ;
}

void time_integrator_free(time_integrator_t* integ)
{
  // Kill the preconditioner stuff.
  Destroy_SuperNode_Matrix(&integ->precond_L);
  Destroy_CompCol_Matrix(&integ->precond_U);
  StatFree(&integ->precond_stat);
  SUPERLU_FREE(integ->precond_cperm);
  SUPERLU_FREE(integ->precond_rperm);
  supermatrix_factory_free(integ->precond_factory);
  supermatrix_free(integ->precond_mat);
  supermatrix_free(integ->precond_rhs);

  // Kill the CVode stuff.
  N_VDestroy(integ->x);
  CVodeFree(&integ->cvode);

  // Kill the rest.
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  free(integ->name);
  free(integ);
}

char* time_integrator_name(time_integrator_t* integ)
{
  return integ->name;
}

void* time_integrator_context(time_integrator_t* integ)
{
  return integ->context;
}

int time_integrator_order(time_integrator_t* integ)
{
  return integ->order;
}

void time_integrator_step(time_integrator_t* integ, real_t t1, real_t t2, real_t* X)
{
  ASSERT(t2 > t1);
  if (integ->current_time != t1)
  {
    CVodeReInit(integ->cvode, t1, integ->x);
    integ->current_time = t1;
  }
  real_t t;
  CVode(integ->cvode, t2, integ->x, &t, CV_NORMAL);
  ASSERT(fabs(t - t2) < 1e-14);
  integ->current_time = t;
}

