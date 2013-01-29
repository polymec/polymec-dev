#include "core/diffusion_solver.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_parcsr_ls.h"

#ifdef __cplusplus
extern "C" {
#endif

struct diffusion_solver_t
{
  char* name;

  // Linear system and solver.
  HYPRE_IJMatrix A;
  HYPRE_IJVector x, b;
  HYPRE_Solver solver;

  // Work vectors.
  HYPRE_IJVector* work;
  int num_work_vectors;

  // Flag is set to true if the linear system above is initialized.
  bool initialized;

  // Table and array for computing A and b components.
  double_table_t* Aij;
  double* bi;

  // Index space.
  index_space_t* index_space;

  // Context pointer, virtual table.
  void* context;
  diffusion_solver_vtable vtable;
};

static void initialize(diffusion_solver_t* solver)
{
  if (!solver->initialized)
  {
    MPI_Comm comm = solver->index_space->comm;
    int low = solver->index_space->low;
    int high = solver->index_space->high;
    HYPRE_IJMatrixCreate(comm, low, high, low, high, &solver->A);
    HYPRE_IJVectorCreate(comm, low, high, &solver->x);
    HYPRE_IJVectorCreate(comm, low, high, &solver->b);
    HYPRE_ParCSRGMRESCreate(comm, &solver->solver);
    solver->initialized = true;
  }
}

static void create_work_vectors(diffusion_solver_t* solver, int num_vectors)
{
  if (solver->num_work_vectors < num_vectors)
  {
    solver->work = realloc(solver->work, sizeof(HYPRE_IJVector)*num_vectors);
    for (int i = solver->num_work_vectors; i < num_vectors; ++i)
      HYPRE_IJVectorCreate(comm, low, high, &solver->work[i]);
    solver->num_work_vectors = num_vectors;
  }
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

  solver->Aij = double_table_new();
  solver->bi = malloc(sizeof(double)*(index_space->high - index_space->low));

  solver->work = NULL;
  solver->num_work_vectors = 0;

  solver->initialized = false;

  return solver;
}

void diffusion_solver_free(diffusion_solver_t* solver)
{
  if (solver->initialized)
  {
    HYPRE_ParCSRGMRESDestroy(solver->solver);
    HYPRE_IJMatrixDestroy(solver->A);
    HYPRE_IJVectorDestroy(solver->x);
    HYPRE_IJVectorDestroy(solver->b);
  }

  // Destroy any work vectors.
  for (int i = 0; i < solver->num_work_vectors; ++i)
    HYPRE_IJVectorDestroy(solver->work[i]);
  free(solver->work);

  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);

  double_table_free(solver->Aij);
  free(solver->bi);

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

static inline void copy_array_to_vector(index_space_t* is, double* array, HYPRE_IJVector vector)
{
  int N = is->high - is->low;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = is->low + i;
  HYPRE_IJVectorSetValues(vector, N, indices, array);
}

static inline void add_array_to_vector(index_space_t* is, double* array, HYPRE_IJVector vector)
{
  int N = is->high - is->low;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = is->low + i;
  HYPRE_IJVectorAddToValues(vector, N, indices, array);
}

static inline void copy_vector_to_array(index_space_t* is, HYPRE_IJVector vector, double* array)
{
  int N = is->high - is->low;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = is->low + i;
  HYPRE_IJVectorGetValues(vector, N, indices, array);
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

static inline void solve(diffusion_solver_t* solver, double_table_t* Aij, double* bi, double* xi)
{
  // Copy the values from Aij and bi to our linear system.

  HYPRE_ParCSRMatrix mat;
  HYPRE_IJMatrixGetObject(solver->A, (void**)&mat);
  HYPRE_ParVector X, B;
  HYPRE_IJVectorGetObject(solver->x, (void**)&X);
  HYPRE_IJVectorGetObject(solver->b, (void**)&B);
  HYPRE_GMRESSetup(solver->solver, (HYPRE_Matrix)mat, (HYPRE_Vector)B, (HYPRE_Vector)X);
  HYPRE_GMRESSolve(solver->solver, (HYPRE_Matrix)mat, (HYPRE_Vector)B, (HYPRE_Vector)X);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, xi);
}

void diffusion_solver_euler(diffusion_solver_t* solver,
                            double t1, double* sol1,
                            double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;
  int N = solver->index_space->high - solver->index_space->low;
  double si[N];

  // Make sure the solver is initialized.
  initialize(solver);

  // A -> diffusion matrix at time t2.
  compute_diff_matrix(solver, solver->Aij, t2);

  // Compute the source at time t2.
  compute_source_vector(solver, si, t2);

  // Apply boundary conditions to the system.
  for (int i = 0; i < N; ++i)
    solver->bi[i] = 0.0;
  apply_bcs(solver, solver->Aij, solver->bi, t2);

  // A -> I - dt*A.
  double_table_val_pos pos = double_table_start(solver->Aij);
  int i, j;
  double Aij;
  while (double_table_next(solver->Aij, &pos, &i, &j, &Aij))
  {
    if (i == j)
      double_table_insert(solver->Aij, i, j, 1.0 - dt * Aij);
    else
      double_table_insert(solver->Aij, i, j, -dt * Aij);
  }

  // What we've done to the LHS we must do to the RHS.
  for (int i = 0; i < N; ++i)
    solver->bi[i] *= -dt;

  // b += sol1 + dt * source.
  for (int i = 0; i < N; ++i)
    solver->bi[i] += sol1[i] + dt * si[i];

  // Solve the linear system.
  solve(solver, solver->Aij, solver->bi, sol2);
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

  // In addition to b and x, we need 3 work vectors, all described below.
  create_work_vectors(solver, 3);

  // Parameters for the TGA algorithm.
  double a = 2.0 - sqrt(2.0);
//  double b = a - 0.5;
  double r1 = (2*a - 1.0) / (a + sqrt(a*a - 4*a + 2.0));
  double r2 = (2*a - 1.0) / (a - sqrt(a*a - 4*a + 2.0));

  // A -> diffusion matrix at time t2.
  compute_diff_matrix(solver, solver->Aij, t2);

  // Apply boundary conditions to the system.
  for (int i = 0; i < N; ++i)
    solver->bi[i] = 0.0;
  apply_bcs(solver, solver->Aij, solver->bi, t2);

  // A note about the "b" vector above: It is intended to ensure that 
  // a solution to the equation A*x = b satisfies boundary conditions.
  // However, we use it in this function as though the diffusion operator 
  // applied to a solution D(x) = A*x + b. This means that we have to 
  // flip its sign.
  for (int i = 0; i < N; ++i)
    solver->bi[i] *= -1.0;

  //-------------------------------------------
  // Construct e, the RHS for the first solve.
  //-------------------------------------------

  // b -> (1.0 - a) * dt * b.
  for (int i = 0; i < N; ++i)
    solver->bi[i] *= (1.0 - a) * dt;

  // e = [I + (1-a) * dt * A] * sol1.
  double ei[N];
  for (int i = 0; i < N; ++i)
    ei[N] = 0.0;
  double_table_val_pos pos = double_table_start(solver->Aij);
  int i, j;
  double Aij;
  while (double_table_next(solver->Aij, &pos, &i, &j, &Aij))
  {
    ei[i] += (1.0 - a) * dt * Aij * sol1[j];
    if (i == j)
      ei[i] += sol1[j];
  }

  // Compute the source at t1.
  double si[N];
  compute_source_vector(solver, si, t1);

  // e -> e + 0.5 * dt * source(t1).
  for (int i = 0; i < N; ++i)
    ei[i] += 0.5 * dt * si[i];

  // Now compute the source at t2. 
  compute_source_vector(solver, si, t2);

  // Transform A  ->  I - 2 * (a - 0.5) * dt * A.
  HYPRE_IJMatrixShift(solver->A, -1.0);
  HYPRE_IJMatrixScale(solver->A, -2.0 * (a - 0.5) / (1.0 - a));
  HYPRE_IJMatrixShift(solver->A, 1.0);
  HYPRE_IJVectorScale(solver->b, -2.0 * (a - 0.5) / (1.0 - a));

  // Compute [I - 2.0 * (a - 0.5) * dt * A] * source(t2), and 
  // store it in work vector 2.
  HYPRE_IJMatrixMultAdd(solver->A, solver->work[0], solver->b, solver->work[2]);

  // e -> e + 0.5 * dt * [I - 2.0 * (a - 0.5) * dt * A * source(t2)
  // (result stored in work vector 1).
  HYPRE_IJVectorAXPY(solver->work[1], 0.5 * dt, solver->work[2]);

  //--------------------------------------------------------
  // Now we have computed e and stored it in work vector 1.
  //--------------------------------------------------------

  // Next, transform (I - 2 * (a - 0.5) * dt * A) to (I - r2 * dt * A).
  HYPRE_IJMatrixShift(solver->A, -1.0);
  HYPRE_IJMatrixScale(solver->A, r2 / (2.0 * (a - 0.5)));
  HYPRE_IJMatrixShift(solver->A, 1.0);
  HYPRE_IJVectorScale(solver->b, r2 / (2.0 * (a - 0.5)));

  // Do the first solve: (I - r2 * dt * A) * v = e. Recall that b contains the 
  // boundary condition information for A, so it needs to be moved to the 
  // right hand side.
  HYPRE_IJVectorAXPY(solver->work[1], -1.0, solver->b); // Move b to RHS
  solve(solver, solver->A, solver->work[1], solver->work[0]); // Solve!
  // The solution v is now stored in work vector 0.

  // Now transform (I - r2 * dt * A) -> (I - r1 * dt * A).
  HYPRE_IJMatrixShift(solver->A, -1.0);
  HYPRE_IJMatrixScale(solver->A, r1 / r2);
  HYPRE_IJMatrixShift(solver->A, 1.0);
  HYPRE_IJVectorScale(solver->b, r1 / r2);

  // Do the second solve. As in the first, we need to move b over to the RHS.
  HYPRE_IJVectorAXPY(solver->work[0], -1.0, solver->b); // Move b to RHS
  solve(solver, solver->A, solver->work[0], solver->x);
  // Solution is now stored in x.

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, sol2);
}

#ifdef __cplusplus
}
#endif

