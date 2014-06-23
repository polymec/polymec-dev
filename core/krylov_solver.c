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

#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"
#include "sundials/sundials_sptfqmr.h"

#include "core/krylov_solver.h"
#include "core/sundials_helpers.h"

typedef enum
{
  GMRES, BICGSTAB, TFQMR
} solver_type;

struct krylov_solver_t 
{
  solver_type type;
  MPI_Comm comm;
  char* name;
  void* context;
  krylov_solver_vtable vtable;
  real_t delta;
  int N, max_krylov_dim, max_restarts;

  SpgmrMem gmres;
  SpbcgMem bicgstab;
  SptfqmrMem tfqmr;

  preconditioner_t* precond;

  N_Vector X, B;
};

krylov_solver_t* gmres_krylov_solver_new(MPI_Comm comm, 
                                         void* context,
                                         krylov_solver_vtable vtable,
                                         int N,
                                         int max_krylov_dim,
                                         int max_restarts)
{
  ASSERT(vtable.ax != NULL);
  ASSERT(N > 0);
  ASSERT(max_krylov_dim >= 3);
  ASSERT(max_restarts >= 0);

  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup("GMRES");
  solver->type = GMRES;
  solver->comm = comm;
  solver->context = context;
  solver->vtable = vtable;
  solver->N = N;
  solver->max_krylov_dim = max_krylov_dim;
  solver->max_restarts = max_restarts;
  solver->delta = 1e-8;

  solver->X = N_VNew(comm, N);
  solver->B = N_VNew(comm, N);
  solver->gmres = SpgmrMalloc(max_krylov_dim, solver->X);
  solver->bicgstab = NULL;
  solver->tfqmr = NULL;
  solver->precond = NULL;

  return solver;
}
  
krylov_solver_t* bicgstab_krylov_solver_new(MPI_Comm comm,
                                            void* context,
                                            krylov_solver_vtable vtable,
                                            int N,
                                            int max_krylov_dim)
{
  ASSERT(vtable.ax != NULL);
  ASSERT(N > 0);
  ASSERT(max_krylov_dim >= 3);

  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup("BiCGSTAB");
  solver->type = BICGSTAB;
  solver->comm = comm;
  solver->context = context;
  solver->vtable = vtable;
  solver->N = N;
  solver->max_krylov_dim = max_krylov_dim;
  solver->delta = 1e-8;

  solver->X = N_VNew(comm, N);
  solver->B = N_VNew(comm, N);
  solver->gmres = NULL;
  solver->bicgstab = SpbcgMalloc(max_krylov_dim, solver->X);
  solver->tfqmr = NULL;
  solver->precond = NULL;

  return solver;
}
  
krylov_solver_t* tfqmr_krylov_solver_new(MPI_Comm comm,
                                         void* context,
                                         krylov_solver_vtable vtable,
                                         int N,
                                         int max_krylov_dim)
{
  ASSERT(vtable.ax != NULL);
  ASSERT(N > 0);
  ASSERT(max_krylov_dim >= 3);

  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup("TFQMR");
  solver->type = TFQMR;
  solver->comm = comm;
  solver->context = context;
  solver->vtable = vtable;
  solver->N = N;
  solver->max_krylov_dim = max_krylov_dim;
  solver->delta = 1e-8;

  solver->X = N_VNew(comm, N);
  solver->B = N_VNew(comm, N);
  solver->gmres = NULL;
  solver->bicgstab = NULL;
  solver->tfqmr = SptfqmrMalloc(max_krylov_dim, solver->X);
  solver->precond = NULL;

  return solver;
}
 
void krylov_solver_free(krylov_solver_t* solver)
{
  if (solver->type == GMRES)
  {
    SpgmrFree(solver->gmres);
  }
  else if (solver->type == BICGSTAB)
  {
    SpbcgFree(solver->bicgstab);
  }
  else
  {
    SptfqmrFree(solver->tfqmr);
  }
  N_VDestroy(solver->X);
  N_VDestroy(solver->B);
  polymec_free(solver->name);
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  polymec_free(solver);
}

void* krylov_solver_context(krylov_solver_t* solver)
{
  return solver->context;
}

void krylov_solver_set_tolerance(krylov_solver_t* solver, real_t delta)
{
  ASSERT(delta > 0.0);
  solver->delta = delta;
}

void krylov_solver_set_preconditioner(krylov_solver_t* solver, preconditioner_t* precond)
{
  solver->precond = precond;
}

// This implements the A*X function in terms understandable to Sundials.
static int krylov_ax(void* solver_ptr, N_Vector x, N_Vector Ax)
{
  krylov_solver_t* solver = solver_ptr;
  return solver->vtable.ax(solver->context, NV_DATA(x), NV_DATA(Ax), solver->N);
}

real_t* krylov_solver_vector(krylov_solver_t* solver)
{
  real_t* X = polymec_malloc(sizeof(real_t) * solver->N);
  memset(X, 0, sizeof(real_t) * solver->N);
  return X;
}

// This is a proxy to the preconditioner solve function.
static int pc_solve(void* context, N_Vector r, N_Vector z, int pc_type)
{
  if (pc_type != PREC_NONE)
  {
    krylov_solver_t* solver = context;
    memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * NV_LOCLENGTH(r));
    bool result = preconditioner_solve(solver->precond, NV_DATA(z));
    return (result) ? 0 : -1;
  }
  return 0;
}

bool krylov_solver_solve(krylov_solver_t* solver, real_t* X, real_t* res_norm, 
                         int* num_iters, int* num_precond)
{
  // Copy data into the X vector.
  memcpy(NV_DATA(solver->X), X, sizeof(real_t) * solver->N);

  // Compute the right hand.
  if (solver->vtable.b != NULL)
    solver->vtable.b(solver->context, NV_DATA(solver->B), solver->N);

  int prec_type = (solver->precond != NULL) ? PREC_RIGHT : PREC_NONE;

  bool result;
  if (solver->type == GMRES)
  {
    int stat = SpgmrSolve(solver->gmres, solver, solver->X, 
                          solver->B, prec_type, MODIFIED_GS, solver->delta, 
                          solver->max_restarts, solver, NULL, NULL, krylov_ax, pc_solve, res_norm,
                          num_iters, num_precond);
    result = (stat == SPGMR_SUCCESS);
  }
  else if (solver->type == BICGSTAB)
  {
    int stat = SpbcgSolve(solver->bicgstab, solver, solver->X, 
                          solver->B, prec_type, solver->delta, solver, 
                          NULL, NULL, krylov_ax, pc_solve, res_norm,
                          num_iters, num_precond);
    result = (stat == SPBCG_SUCCESS);
  }
  else
  {
    int stat = SptfqmrSolve(solver->tfqmr, solver, solver->X, 
                            solver->B, prec_type, solver->delta, solver, 
                            NULL, NULL, krylov_ax, pc_solve, res_norm,
                            num_iters, num_precond);
    result = (stat == SPTFQMR_SUCCESS);
  }
  
  // Copy the data out.
  memcpy(X, NV_DATA(solver->X), sizeof(real_t) * solver->N);
  return result;
}

