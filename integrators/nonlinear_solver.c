// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
  int block_size;
  void (*dtor)(void*);
  mf_nonlinear_solver_type_t solver_type;
  adj_graph_t* graph;
  void* kinsol;

  // Preconditioning stuff.
  SuperMatrix precond_mat, precond_rhs, precond_L, precond_U;
  int *precond_rperm, *precond_cperm;
};

nonlinear_solver_t* nonlinear_solver_new(const char* name, 
                                         void* context,
                                         KINSysFn F,
                                         void (*dtor)(void*),
                                         adj_graph_t* graph,
                                         mf_nonlinear_solver_type_t type)
{
  ASSERT(F != NULL);
  nonlinear_solver_t* solver = malloc(sizeof(nonlinear_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->graph = graph;
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
  free_preconditioner(solver);
  KINFree(&solver->kinsol);
  if ((solver->dtor != NULL) && (solver->context != NULL))
    solver->dtor(solver->context);
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
                                  
