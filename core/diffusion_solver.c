#include "core/diffusion_solver.h"
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

  // Index range.
  int ilow, ihigh;

  // Context pointer, virtual table.
  void* context;
  diffusion_solver_vtable vtable;
};

static void initialize(diffusion_solver_t* solver)
{
  if (!solver->initialized)
  {
    solver->vtable.create_matrix(solver->context, &solver->A);
    solver->vtable.create_vector(solver->context, &solver->x);
    solver->vtable.create_vector(solver->context, &solver->b);
    solver->vtable.create_solver(solver->context, &solver->solver);
    HYPRE_IJVectorGetLocalRange(solver->x, &solver->ilow, &solver->ihigh);
    solver->initialized = true;
  }
}

static void create_work_vectors(diffusion_solver_t* solver, int num_vectors)
{
  if (solver->num_work_vectors < num_vectors)
  {
    solver->work = realloc(solver->work, sizeof(HYPRE_IJVector)*num_vectors);
    for (int i = solver->num_work_vectors; i < num_vectors; ++i)
      solver->vtable.create_vector(solver->context, &solver->work[i]);
    solver->num_work_vectors = num_vectors;
  }
}

diffusion_solver_t* diffusion_solver_new(const char* name, 
                                         void* context,
                                         diffusion_solver_vtable vtable)
{
  ASSERT(name != NULL);
  ASSERT(vtable.create_matrix != NULL);
  ASSERT(vtable.create_vector != NULL);
  ASSERT(vtable.create_solver != NULL);
  ASSERT(vtable.compute_diffusion_matrix != NULL);
  ASSERT(vtable.apply_bcs != NULL);

  diffusion_solver_t* solver = malloc(sizeof(diffusion_solver_t));
  solver->name = strdup(name);
  solver->context = context;
  solver->vtable = vtable;

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

static inline void copy_array_to_vector(double* array, HYPRE_IJVector vector)
{
  int ilow, ihigh;
  HYPRE_IJVectorGetLocalRange(vector, &ilow, &ihigh);
  int N = ihigh - ilow;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = ilow + i;
  HYPRE_IJVectorSetValues(vector, N, indices, array);
}

static inline void add_array_to_vector(double* array, HYPRE_IJVector vector)
{
  int ilow, ihigh;
  HYPRE_IJVectorGetLocalRange(vector, &ilow, &ihigh);
  int N = ihigh - ilow;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = ilow + i;
  HYPRE_IJVectorAddToValues(vector, N, indices, array);
}

static inline void copy_vector_to_array(HYPRE_IJVector vector, double* array)
{
  int ilow, ihigh;
  HYPRE_IJVectorGetLocalRange(vector, &ilow, &ihigh);
  int N = ihigh - ilow;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = ilow + i;
  HYPRE_IJVectorGetValues(vector, N, indices, array);
}

static inline void compute_diff_matrix(diffusion_solver_t* solver, HYPRE_IJMatrix A, double t)
{
  solver->vtable.compute_diffusion_matrix(solver->context, A, t);
}

static inline void compute_source_vector(diffusion_solver_t* solver, HYPRE_IJVector source, double t)
{
  solver->vtable.compute_source_vector(solver->context, source, t);
}

static inline void apply_bcs(diffusion_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, double t)
{
  solver->vtable.apply_bcs(solver->context, A, b, t);
}

static inline void solve(diffusion_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_ParCSRMatrix mat;
  HYPRE_IJMatrixGetObject(A, (void**)&mat);
  HYPRE_ParVector X, B;
  HYPRE_IJVectorGetObject(x, (void**)&X);
  HYPRE_IJVectorGetObject(b, (void**)&B);
  HYPRE_GMRESSetup(solver->solver, (HYPRE_Matrix)mat, (HYPRE_Vector)B, (HYPRE_Vector)X);
  HYPRE_GMRESSolve(solver->solver, (HYPRE_Matrix)mat, (HYPRE_Vector)B, (HYPRE_Vector)X);
}

void diffusion_solver_euler(diffusion_solver_t* solver,
                            double t1, double* sol1,
                            double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;

  // Make sure the solver is initialized.
  initialize(solver);

  // A -> diffusion matrix at time t2.
  compute_diff_matrix(solver, solver->A, t2);

  // Compute the source at time t2.
  compute_source_vector(solver, solver->x, t2);

  // Apply boundary conditions to the system.
  int N = solver->ihigh - solver->ilow;
  int indices[N];
  double values[N];
  for (int i = 0; i < N; ++i)
  {
    indices[i] = solver->ilow + i;
    values[N] = 0.0;
  }
  HYPRE_IJVector_SetValues(solver->b, N, indices, values);
  apply_bcs(solver, solver->A, solver->b, t2);

  // A -> I - dt*A.
  HYPRE_IJMatrixScale(solver->A, -dt);
  HYPRE_IJMatrixShift(solver->A, 1.0);

  // What we've done to the LHS we must do to the RHS.
  HYPRE_IJVectorScale(solver->b, -dt);

  // b += sol1 + dt * source.
  add_array_to_vector(sol1, solver->b);
  HYPRE_IJVectorAXPY(solver->b, dt, solver->x);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, sol2);
}

void diffusion_solver_tga(diffusion_solver_t* solver,
                          double t1, double* sol1, 
                          double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;

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
  compute_diff_matrix(solver, solver->A, t2);

  // Apply boundary conditions to the system.
  HYPRE_IJVectorSet(solver->b, 0.0);
  apply_bcs(solver, solver->A, solver->b, t2);

  // A note about the "b" vector above: It is intended to ensure that 
  // a solution to the equation A*x = b satisfies boundary conditions.
  // However, we use it in this function as though the diffusion operator 
  // applied to a solution D(x) = A*x + b. This means that we have to 
  // flip its sign.
  HYPRE_IJVectorScale(solver->b, -1.0);

  //-------------------------------------------
  // Construct e, the RHS for the first solve.
  //-------------------------------------------

  // A -> (1-a) * dt * A.
  HYPRE_IJMatrixScale(solver->A, (1.0 - a) * dt);
  HYPRE_IJVectorScale(solver->b, (1.0 - a) * dt);
  // A -> A + I.
  HYPRE_IJMatrixShift(solver->A, 1.0);
  // e = [I + (1-a) * dt * A] * sol1 (stored in work vector 1).
  copy_array_to_vector(sol1, solver->work[0]);
  HYPRE_IJMatrixMultAdd(solver->A, solver->work[0], solver->b, solver->work[1]);

  // Compute the source at t1. We'll store the result in work vector 0.
  compute_source_vector(solver, solver->work[0], t1);

  // e -> e + 0.5 * dt * source(t1), stored in work vector 1.
  HYPRE_IJVectorAXPY(solver->work[1], 0.5 * dt, solver->work[0]);

  // Now compute the source at t2. We'll store the result in work vector 0.
  compute_source_vector(solver, solver->work[0], t2);

  // Transform I + (1 - a) * dt * A  ->  I - 2 * (a - 0.5) * dt * A.
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

