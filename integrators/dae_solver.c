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
#include "integrators/dae_solver.h"

struct dae_solver_t 
{
  void* context;
  char* name;
  dae_solver_vtable vtable;
  int N;
};

dae_solver_t* dae_solver_new(const char* name, 
                             void* context,
                             dae_solver_vtable vtable, 
                             int N)
{
  ASSERT(vtable.eval_dae != NULL);
  ASSERT(N > 0);

  dae_solver_t* solver = malloc(sizeof(dae_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->N = N;

  return solver;
}

void dae_solver_free(dae_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  free(solver->name);
  free(solver);
}

char* dae_solver_name(dae_solver_t* solver)
{
  return solver->name;
}

void* dae_solver_context(dae_solver_t* solver)
{
  return solver->context;
}

void dae_solver_solve(dae_solver_t* solver,
                      double t,
                      double* X)
{
  ASSERT(X != NULL);
  solver->vtable.solve(solver->context, t, solver->N, X);
}
                                  
