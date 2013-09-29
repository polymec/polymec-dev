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

//#include "core/hypre_helpers.h"
#include "integrators/elliptic_solver.h"

struct elliptic_solver_t
{
  char* name;

  // Linear system and solver.
  void *A, *x, *b, *solver, *pc;
#if 0
  HYPRE_IJMatrix A;
  HYPRE_IJVector x, b;
  HYPRE_Solver solver;
//  HYPRE_Solver pc; // Preconditioner
#endif

  bool matrix_created;

  // Index space.
  index_space_t* index_space;

  // Context pointer, virtual table.
  void* context;
  elliptic_solver_vtable vtable;
};

static void initialize(elliptic_solver_t* solver)
{
#if 0
  solver->A = HYPRE_IJMatrixNew(solver->index_space);
  solver->x = HYPRE_IJVectorNew(solver->index_space);
  solver->b = HYPRE_IJVectorNew(solver->index_space);
  HYPRE_ParCSRHybridCreate(&solver->solver);
  HYPRE_ParCSRHybridSetSolverType(solver->solver, 2);
  HYPRE_ParCSRHybridSetKDim(solver->solver, 5);

  // Set up the preconditioner.
//  HYPRE_ParaSailsCreate(solver->index_space->comm, &solver->pc);
//  HYPRE_ParaSailsSetSym(solver->pc, 0);
//  HYPRE_GMRESSetPrecond(solver->solver, 
//                        (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
//                        (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup,
//                        solver->pc);
#endif
}

elliptic_solver_t* elliptic_solver_new(const char* name, 
                                       void* context,
                                       elliptic_solver_vtable vtable,
                                       index_space_t* index_space)
{
  ASSERT(name != NULL);
  ASSERT(vtable.compute_operator_matrix != NULL);
  ASSERT(vtable.apply_bcs != NULL);
  ASSERT(index_space != NULL);

  elliptic_solver_t* solver = malloc(sizeof(elliptic_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->index_space = index_space;
  solver->matrix_created = false;

  solver->A = NULL;
  solver->b = NULL;
  solver->x = NULL;
//  solver->pc = NULL;
  solver->solver = NULL;

  // Make sure the solver is initialized.
  initialize(solver);

  return solver;
}

void elliptic_solver_free(elliptic_solver_t* solver)
{
#if 0
//  if (solver->pc != NULL)
//    HYPRE_ParaSailsDestroy(solver->pc);
  HYPRE_ParCSRHybridDestroy(solver->solver);
  HYPRE_IJMatrixDestroy(solver->A);
  HYPRE_IJVectorDestroy(solver->x);
  HYPRE_IJVectorDestroy(solver->b);
#endif

  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);

  free(solver->name);
  free(solver);
}

char* elliptic_solver_name(elliptic_solver_t* solver)
{
  return solver->name;
}

void* elliptic_solver_context(elliptic_solver_t* solver)
{
  return solver->context;
}

#if 0
static inline void copy_vector_to_array(index_space_t* is, HYPRE_IJVector vector, double* array)
{
  int N = is->high - is->low;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = is->low + i;
  HYPRE_IJVectorGetValues(vector, N, indices, array);
}
#endif

#if 0
static inline void solve(elliptic_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
//HYPRE_IJMatrixPrint(solver->A, "A.txt");
  HYPRE_ParCSRMatrix Aobj;
  int err = HYPRE_IJMatrixGetObject(solver->A, (void**)&Aobj);
  ASSERT(err == 0);
  ASSERT(Aobj != NULL);
  HYPRE_ParVector xobj, bobj;
  err = HYPRE_IJVectorGetObject(solver->x, (void**)&xobj);
  ASSERT(err == 0);
  ASSERT(xobj != NULL);
  err = HYPRE_IJVectorGetObject(solver->b, (void**)&bobj);
  ASSERT(err == 0);
  ASSERT(bobj != NULL);

  err = HYPRE_ParCSRHybridSetup(solver->solver, Aobj, bobj, xobj);
  ASSERT(err == 0);
  err = HYPRE_ParCSRHybridSolve(solver->solver, Aobj, bobj, xobj);
  if (err == HYPRE_ERROR_CONV)
    polymec_error("elliptic_solver: Solve did not converge.");

  HYPRE_ClearAllErrors();
}
#endif

void elliptic_solver_solve(elliptic_solver_t* solver,
                           double t, double* solution)
{
  int N = solver->index_space->high - solver->index_space->low;

  // A -> operator matrix at time t.
  double_table_t* A = double_table_new();
  solver->vtable.compute_operator_matrix(solver, A, t);

  // Compute the source at time t.
  double b[N];
  solver->vtable.compute_source_vector(solver, b, t);

  // Apply boundary conditions to the system.
  solver->vtable.apply_bcs(solver, A, b, t);

#if 0
  // Set up the linear system.
  if (!solver->matrix_created)
  {
    HYPRE_IJMatrixSetRowSizesFromTable(solver->A, solver->index_space, A);
    solver->matrix_created = true;
  }
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, A);
  HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, b);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->index_space, solver->x, solution);
#endif

  // Clean up.
  double_table_free(A);
}

