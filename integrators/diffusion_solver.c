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
#include "core/index_space.h"
#include "core/linear_algebra.h"
#include "integrators/diffusion_solver.h"

struct diffusion_solver_t
{
  char* name;

  // Linear system and solver.
  void *A, *x, *b, *solver;
#if 0
  HYPRE_IJMatrix A;
  HYPRE_IJVector x, b;
  HYPRE_Solver solver;
#endif

  // Flag is set to true if the linear system above is initialized.
  bool initialized;

  // Index space.
  index_space_t* index_space;

  // Context pointer, virtual table.
  void* context;
  diffusion_solver_vtable vtable;
};

static void initialize(diffusion_solver_t* solver)
{
#if 0
  if (!solver->initialized)
  {
    solver->A = HYPRE_IJMatrixNew(solver->index_space);
    solver->x = HYPRE_IJVectorNew(solver->index_space);
    solver->b = HYPRE_IJVectorNew(solver->index_space);
    HYPRE_ParCSRHybridCreate(&solver->solver);
    HYPRE_ParCSRHybridSetSolverType(solver->solver, 2);
    HYPRE_ParCSRHybridSetKDim(solver->solver, 5);
    solver->initialized = true;
  }
#endif
}

diffusion_solver_t* diffusion_solver_new(const char* name, 
                                         void* context,
                                         diffusion_solver_vtable vtable,
                                         index_space_t* index_space)
{
  ASSERT(name != NULL);
  ASSERT(vtable.compute_diffusion_matrix != NULL);
  ASSERT(vtable.apply_bcs != NULL);
  ASSERT(index_space != NULL);

  diffusion_solver_t* solver = malloc(sizeof(diffusion_solver_t));
  solver->name = strdup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->index_space = index_space;

  solver->initialized = false;

  return solver;
}

void diffusion_solver_free(diffusion_solver_t* solver)
{
#if 0
  if (solver->initialized)
  {
    HYPRE_ParCSRHybridDestroy(solver->solver);
    HYPRE_IJMatrixDestroy(solver->A);
    HYPRE_IJVectorDestroy(solver->x);
    HYPRE_IJVectorDestroy(solver->b);
  }
#endif

  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);

  free(solver->name);
  free(solver);
}

char* diffusion_solver_name(diffusion_solver_t* solver)
{
  return solver->name;
}

void* diffusion_solver_context(diffusion_solver_t* solver)
{
  return solver->context;
}

static inline void compute_diff_matrix(diffusion_solver_t* solver, double_table_t* A, double t)
{
  solver->vtable.compute_diffusion_matrix(solver->context, A, t);
}

static inline void compute_source_vector(diffusion_solver_t* solver, double* source, double t)
{
  solver->vtable.compute_source_vector(solver->context, source, t);
}

static inline void apply_bcs(diffusion_solver_t* solver, double_table_t* A, double* b, double t)
{
  solver->vtable.apply_bcs(solver->context, A, b, t);
}

#if 0
static inline void solve(diffusion_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
//HYPRE_IJMatrixPrint(solver->A, "A");
//HYPRE_IJVectorPrint(solver->b, "b");
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
    polymec_error("diffusion_solver: Solve did not converge.");

  HYPRE_ClearAllErrors();
}
#endif

void diffusion_solver_euler(diffusion_solver_t* solver,
                            double t1, double* sol1,
                            double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;
  int N = solver->index_space->high - solver->index_space->low;

  // Make sure the solver is initialized.
  initialize(solver);

  // A -> diffusion matrix at time t2.
  double_table_t* A = double_table_new();
  compute_diff_matrix(solver, A, t2);

  // Apply boundary conditions to the system.
  double b[N];
  for (int i = 0; i < N; ++i)
    b[i] = 0.0;
  apply_bcs(solver, A, b, t2);

  // L <- I - dt*A.
  double_table_cell_pos_t pos = double_table_start(A);
  int i = 0, j = 0;
  double Aij;
  double_table_t* L = double_table_new();
  while (double_table_next_cell(A, &pos, &i, &j, &Aij))
  {
//printf("A(%d, %d) = %g\n", i, j, Aij);
    if (i == j)
      double_table_insert(L, i, j, 1.0 - dt * Aij);
    else
      double_table_insert(L, i, j, -dt * Aij);
  }

  // Compute the source at time t2.
  double si[N];
  compute_source_vector(solver, si, t2);

  // b = -dt*(bc terms) + sol1 + dt * source.
  for (int i = 0; i < N; ++i)
    b[i] = -dt * b[i] + sol1[i] + dt * si[i];

#if 0
  // Set up the linear system.
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, L);
//HYPRE_IJMatrixPrint(solver->A, "L");
  HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, b);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, sol2);
//printf("x = ");
//vector_fprintf(sol2, N, stdout);
//printf("\n");
#endif

  // Clean up.
  double_table_free(L);
  double_table_free(A);
}

void diffusion_solver_tga(diffusion_solver_t* solver,
                          double t1, double* sol1, 
                          double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;
  int N = solver->index_space->high - solver->index_space->low;

  // Make sure the solver is initialized.
  initialize(solver);

  // Parameters for the TGA algorithm.
  double a = 2.0 - sqrt(2.0);
//  double b = a - 0.5;
  double r1 = (2*a - 1.0) / (a + sqrt(a*a - 4*a + 2.0));
  double r2 = (2*a - 1.0) / (a - sqrt(a*a - 4*a + 2.0));

  // A -> Components of diffusion matrix at time t2.
  double_table_t* A = double_table_new();
  compute_diff_matrix(solver, A, t2);

  // Apply boundary conditions to the system.
  double b[N];
  for (int i = 0; i < N; ++i)
    b[i] = 0.0;
  apply_bcs(solver, A, b, t2);

  //-------------------------------------------
  // Construct e, the RHS for the first solve.
  //-------------------------------------------

  // Compute the source at t1 and t2.
  double s1[N], s2[N];
  compute_source_vector(solver, s1, t1);
  compute_source_vector(solver, s2, t2);

  // e = [I + (1 - a) * dt * A] * sol1 + 
  //     0.5 * dt * [s1 + [I - 2*(a - 0.5) * dt * A] * s2].
  // Also, M1 = I - r2 * dt * A and M2 = I - r1 * dt * A.
  double e[N];
  for (int i = 0; i < N; ++i)
    e[i] = 0.0;
  double_table_cell_pos_t pos = double_table_start(A);
  int i = 0, j = 0;
  double Aij;
  double_table_t* M1 = double_table_new();
  double_table_t* M2 = double_table_new();
  while (double_table_next_cell(A, &pos, &i, &j, &Aij))
  {
    // Compute local indices for the global indices i and j.
    int ii = i - solver->index_space->low;
    int jj = j - solver->index_space->low;

    // Add the off-diagonal terms for e.
    e[ii] += (1.0 - a) * dt * Aij * sol1[jj] -
             (a - 0.5) * dt * dt * Aij * s2[jj];

    // Set up M1 and M2.
    if (i == j)
    {
      double_table_insert(M1, i, j, 1.0 - r2 * dt * Aij);
      double_table_insert(M2, i, j, 1.0 - r1 * dt * Aij);
    }
    else
    {
      double_table_insert(M1, i, j, -r2 * dt * Aij);
      double_table_insert(M2, i, j, -r1 * dt * Aij);
    }
  }

  // Add the purely diagonal terms to e.
  for (int i = 0; i < N; ++i)
  {
    e[i] += sol1[i] + 0.5 * dt * (s1[i] + s2[i]);

    // Add boundary terms by subtracting the corresponding component of b.
    double bi = b[i];
    e[i] += -(1.0 - a) * dt * bi +  
             (a - 0.5) * dt * dt * bi + 
             -r2 * dt * bi; 
  }
  
#if 0
  // Now solve the linear system M1 * v = e.
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, M1);
  HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, e);
  solve(solver, solver->A, solver->b, solver->x);
#endif
  double v[N];
#if 0
  HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, v);
#endif

  // Add the boundary terms to v.
  pos = double_table_start(A);
  for (int i = 0; i < N; ++i)
    v[i] -= r1 * dt * b[i];

#if 0
  // Now set up the linear system M2 * sol2 = v.
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, M2);
  HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, v);
  solve(solver, solver->A, solver->b, solver->x);
  HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, sol2);
#endif

  // Clean up.
  double_table_free(A);
  double_table_free(M1);
  double_table_free(M2);
}

void diffusion_solver_crank_nicolson(diffusion_solver_t* solver,
                                     double t1, double* sol1,
                                     double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;
  int N = solver->index_space->high - solver->index_space->low;

  // Make sure the solver is initialized.
  initialize(solver);

  // A -> diffusion matrix at time t2.
  double_table_t* A = double_table_new();
  compute_diff_matrix(solver, A, t2);

  // Apply boundary conditions to the system.
  double b[N];
  for (int i = 0; i < N; ++i)
    b[i] = 0.0;
  apply_bcs(solver, A, b, t2);

  // L <- I - 0.5*dt*A, rhs <- (I + 0.5*dt*A)*x
  double_table_cell_pos_t pos = double_table_start(A);
  int i = 0, j = 0;
  double Aij;
  double_table_t* L = double_table_new();
  double rhs[N];
  memset(rhs, 0, sizeof(double)*N);
  while (double_table_next_cell(A, &pos, &i, &j, &Aij))
  {
//printf("A(%d, %d) = %g\n", i, j, Aij);
    if (i == j)
      double_table_insert(L, i, j, 1.0 - 0.5 * dt * Aij);
    else
      double_table_insert(L, i, j, -0.5 * dt * Aij);

    // Add the off-diagonal terms to the RHS.
    int ii = i - solver->index_space->low;
    int jj = j - solver->index_space->low;
    rhs[ii] += 0.5 * dt * Aij * sol1[jj];
  }

  // Compute the source at times t1 and t2.
  double s1[N], s2[N];
  compute_source_vector(solver, s1, t1);
  compute_source_vector(solver, s2, t2);

  // Add the purely diagonal terms to the right hand side.
  for (int i = 0; i < N; ++i)
  {
    rhs[i] += sol1[i] + 0.5 * dt * (s1[i] + s2[i]);

    // Add boundary terms by subtracting the corresponding component of b.
    double bi = b[i];
    rhs[i] += -dt * bi;
  }
 
#if 0
  // Set up the linear system.
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, L);
//HYPRE_IJMatrixPrint(solver->A, "L");
  HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, rhs);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, sol2);
#endif

  // Clean up.
  double_table_free(L);
  double_table_free(A);
}

