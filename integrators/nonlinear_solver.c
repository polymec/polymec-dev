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
#include "integrators/nonlinear_solver.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

// We use a serial version of SuperLU to do preconditioning.
#include "slu_ddefs.h"
#include "supermatrix.h"
#include "slu_util.h"

struct nonlinear_solver_t 
{
  // Parallel stuff.
  int rank, nprocs;

  char* name;
  void* context;
  void* kinsol;
  nonlinear_solver_vtable vtable;
  nonlinear_solver_type_t solver_type;

  // Preconditioning stuff.
  SuperMatrix precond_mat, precond_rhs, precond_L, precond_U;
  int *precond_rperm, *precond_cperm;
};

nonlinear_solver_t* nonlinear_solver_new(const char* name, 
                                         void* context,
                                         nonlinear_solver_vtable vtable,
                                         nonlinear_solver_type_t type)
{
  nonlinear_solver_t* solver = malloc(sizeof(nonlinear_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->solver_type = type;

  return solver;
}

static void free_preconditioner(nonlinear_solver_t* solver)
{
  SUPERLU_FREE(solver->precond_cperm);
  SUPERLU_FREE(solver->precond_rperm);
//  Destroy_CompCol_Matrix(&mf_solver->precond_U);
//  Destroy_SuperNode_Matrix(&mf_solver->precond_L);
  Destroy_SuperMatrix_Store(&solver->precond_rhs);
  Destroy_CompRow_Matrix(&solver->precond_mat);
}

void nonlinear_solver_free(nonlinear_solver_t* solver)
{
//  free_preconditioner(solver);
//  KINFree(&solver->kinsol);
  if ((solver->vtable.dtor != NULL) && (solver->context != NULL))
    solver->vtable.dtor(solver->context);
  free(solver->name);
  free(solver);
}

char* nonlinear_solver_name(nonlinear_solver_t* solver)
{
  return solver->name;
}

void* nonlinear_solver_context(nonlinear_solver_t* solver)
{
  return solver->context;
}

void nonlinear_solver_solve(nonlinear_solver_t* solver,
                            double t,
                            double* X)
{
  ASSERT(X != NULL);
}
                                  
