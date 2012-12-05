#include "core/diffusion_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

struct diffusion_solver_t
{
  char* name;

  // Linear system and solver.
  Mat A;
  Vec x, b;
  KSP ksp;

  // Work vectors.
  Vec* work_vectors;
  int num_work_vectors;

  // Flag is set to true if the linear system above is initialized.
  bool initialized;

  // Context pointer, virtual table.
  void* context;
  diffusion_solver_vtable vtable;
};

static void initialize(diffusion_solver_t* solver)
{
  solver->vtable.create_matrix(solver->context, &solver->A);
  solver->vtable.create_vector(solver->context, &solver->x);
  solver->vtable.create_vector(solver->context, &solver->b);
  solver->vtable.create_ksp(solver->context, &solver->ksp);
}

static void create_work_vectors(diffusion_solver_t* solver, int num_vectors)
{
  if (solver->num_work_vectors < num_vectors)
  {
    solver->work_vectors = realloc(solver->work_vectors, sizeof(Vec)*num_vectors);
    for (int i = solver->num_work_vectors; i < num_vectors; ++i)
      solver->vtable.create_vector(solver->context, &solver->work_vectors[i]);
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
  ASSERT(vtable.create_ksp != NULL);
  ASSERT(vtable.compute_diffusion_matrix != NULL);
  ASSERT(vtable.apply_bcs != NULL);

  diffusion_solver_t* solver = malloc(sizeof(diffusion_solver_t));
  solver->name = strdup(name);
  solver->context = context;
  solver->vtable = vtable;

  solver->work_vectors = NULL;
  solver->num_work_vectors = 0;

  // Make sure the solver is initialized.
  initialize(solver);

  return solver;
}

void diffusion_solver_free(diffusion_solver_t* solver)
{
  KSPDestroy(&solver->ksp);
  MatDestroy(&solver->A);
  VecDestroy(&solver->x);
  VecDestroy(&solver->b);

  // Destroy any work vectors.
  for (int i = 0; i < solver->num_work_vectors; ++i)
    VecDestroy(&solver->work_vectors[i]);

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

static inline void copy_array_to_vector(double* array, Vec vector)
{
  int size;
  VecGetLocalSize(vector, &size);
  double* v;
  VecGetArray(vector, &v);
  memcpy(v, array, sizeof(double)*size);
  VecRestoreArray(vector, &v);
}

static inline void add_array_to_vector(double* array, Vec vector)
{
  int size;
  VecGetLocalSize(vector, &size);
  double* v;
  VecGetArray(vector, &v);
  for (int i = 0; i < size; ++i)
    v[i] += array[i];
  VecRestoreArray(vector, &v);
}

static inline void copy_vector_to_array(Vec vector, double* array)
{
  int size;
  VecGetLocalSize(vector, &size);
  double* v;
  VecGetArray(vector, &v);
  memcpy(array, v, sizeof(double)*size);
  VecRestoreArray(vector, &v);
}

static inline void compute_diff_matrix(diffusion_solver_t* solver, Mat A, double t)
{
  solver->vtable.compute_diffusion_matrix(solver->context, A, t);
}

static inline void compute_source_vector(diffusion_solver_t* solver, Vec source, double t)
{
  solver->vtable.compute_source_vector(solver->context, source, t);
}

static inline void apply_bcs(diffusion_solver_t* solver, Mat A, Vec b, double t)
{
  solver->vtable.apply_bcs(solver->context, A, b, t);
}

static inline void solve(diffusion_solver_t* solver, Mat A, Vec b, Vec x)
{
  KSPSetOperators(solver->ksp, A, A, SAME_NONZERO_PATTERN);
  KSPSolve(solver->ksp, b, x);
}

void diffusion_solver_euler(diffusion_solver_t* solver,
                            double t1, double* sol1,
                            double t2, double* sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;

  // A -> diffusion matrix at time t2.
  compute_diff_matrix(solver, solver->A, t2);

  // Compute the source at time t2.
  compute_source_vector(solver, solver->x, t2);

  // Apply boundary conditions to the system.
  VecSet(solver->b, 0.0);
  apply_bcs(solver, solver->A, solver->b, t2);

  // A -> I - dt*A.
  MatScale(solver->A, -dt);
  MatShift(solver->A, 1.0);

  // What we've done to the LHS we must do to the RHS.
  VecScale(solver->b, -dt);

  // b += sol1 + dt * source.
  add_array_to_vector(sol1, solver->b);
  VecAXPY(solver->b, dt, solver->x);

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
  VecSet(solver->b, 0.0);
  apply_bcs(solver, solver->A, solver->b, t2);

  //-------------------------------------
  // e = [I + (1-a) * dt * A] * sol1.
  //-------------------------------------

  // A -> (1-a) * dt * A.
  MatScale(solver->A, (1.0 - a) * dt);
  VecScale(solver->b, (1.0 - a) * dt);
  // A -> A + I.
  MatShift(solver->A, 1.0);
  // e = [I + (1-a) * dt * A] * sol1 (stored in work vector 1).
  copy_array_to_vector(sol1, solver->work_vectors[0]);
  MatMultAdd(solver->A, solver->work_vectors[0], solver->b, solver->work_vectors[1]);

  // Compute the source at t1. We'll store the result in work vector 0.
  compute_source_vector(solver, solver->work_vectors[0], t1);

  // e -> e + 0.5 * dt * source(t1), stored in work vector 1.
  VecAXPY(solver->work_vectors[1], 0.5 * dt, solver->work_vectors[0]);

  // Now compute the source at t2. We'll store the result in work vector 0.
  compute_source_vector(solver, solver->work_vectors[0], t2);

  // Transform I + (1 - a) * dt * A  ->  I - 2 * (a - 0.5) * dt * A.
  MatShift(solver->A, -1.0);
  MatScale(solver->A, -2.0 * (a - 0.5) / (1.0 - a));
  MatShift(solver->A, 1.0);
  VecScale(solver->b, -2.0 * (a - 0.5) / (1.0 - a));

  // Compute 0.5 * dt * [I - 2.0 * (a - 0.5) * dt * A * source(t2), and 
  // store it in work vector 2.
  MatMult(solver->A, solver->work_vectors[0], solver->work_vectors[2]);

  // e -> e + 0.5 * dt * [I - 2.0 * (a - 0.5) * dt * A * source(t2)
  // (result stored in work vector 1).
  MatMultAdd(solver->A, solver->work_vectors[2], solver->b, solver->work_vectors[1]);

  //--------------------------------------------------------
  // Now we have computed e and stored it in work vector 1.
  //--------------------------------------------------------

  // Next, transform (I - 2 * (a - 0.5) * dt * A) to (I - r2 * dt * A).
  MatShift(solver->A, -1.0);
  MatScale(solver->A, r2 / (2.0 * (a - 0.5)));
  MatShift(solver->A, 1.0);
  VecScale(solver->b, r2 / (2.0 * (a - 0.5)));

  // Do the first solve: (I - r2 * dt * A) v = e. Recall that b contains the 
  // boundary condition information for A, so it needs to be moved to the 
  // right hand side.
  VecAXPY(solver->work_vectors[1], -1.0, solver->b); // Move b to RHS
  solve(solver, solver->A, solver->work_vectors[1], solver->work_vectors[0]); // Solve!
  // The solution is now stored in work vector 0.

  // Now transform (I - r2 * dt * A) -> (I - r1 * dt * A).
  MatShift(solver->A, -1.0);
  MatScale(solver->A, r1 / r2);
  MatShift(solver->A, 1.0);
  VecScale(solver->b, r1 / r2);

  // Do the second solve. As in the first, we need to move b over to the RHS.
  VecAXPY(solver->work_vectors[0], -1.0, solver->b); // Move b to RHS
  solve(solver, solver->A, solver->work_vectors[0], solver->x);
  // Solution is now stored in x.

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, sol2);
}

#ifdef __cplusplus
}
#endif

