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

#include "core/polymec.h"
#include "core/sparse_lin_solver.h"

struct sparse_lin_solver_t
{
  char* name;
  void* context;
  sparse_lin_solver_vtable vtable;

  // Solver metadata.
  double res_norm;
  int nli, nps;

  // Details of last solve outcome.
  char outcome_details[SPARSE_LIN_SOLVER_OUTCOME_DETAILS_MAXLEN];
};

sparse_lin_solver_t* sparse_lin_solver_new(const char* name, void* context,
                                           sparse_lin_solver_vtable vtable)
{
  ASSERT(vtable.solve != NULL);
  sparse_lin_solver_t* solver = malloc(sizeof(sparse_lin_solver_t));
  solver->name = strdup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->res_norm = 0.0;
  solver->nli = 0;
  solver->nps = 0;
  return solver;
}

void sparse_lin_solver_free(sparse_lin_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  free(solver->name);
  free(solver);
}

sparse_lin_solver_outcome_t sparse_lin_solver_solve(sparse_lin_solver_t* solver, void* x, void* b)
{
  solver->res_norm = 0.0;
  solver->nli = 0;
  solver->nps = 0;
  solver->outcome_details[0] = '\0';
  return solver->vtable.solve(solver->context, x, b, &solver->res_norm, &solver->nli, &solver->nps, solver->outcome_details);
}

void sparse_lin_solver_get_info(sparse_lin_solver_t* solver,
                                double* res_l2_norm,
                                int* num_linear_iterations,
                                int* num_precond_solves)
{
  *res_l2_norm = solver->res_norm;
  *num_linear_iterations = solver->nli;
  *num_precond_solves = solver->nps;
}

char* sparse_lin_solver_outcome_details(sparse_lin_solver_t* solver)
{
  return &solver->outcome_details[0];
}

