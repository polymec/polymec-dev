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

#if 0
static inline void copy_array_to_vector(double* array, Vec vector)
{
  int size;
  VecGetLocalSize(vector, &size);
  double* v;
  VecGetArray(vector, &v);
  memcpy(v, array, sizeof(double)*size);
  VecRestoreArray(vector, &v);
}
#endif

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

  // Parameters for the TGA algorithm.
  double a = 2.0 - sqrt(2.0);
  double b = a - 0.5;
  double r1 = (2*a - 1.0) / (a + sqrt(a*a - 4*a + 2.0));
  double r2 = (2*a - 1.0) / (a - sqrt(a*a - 4*a + 2.0));

  // We will do all our work in the new_sol and work vectors.

  // e = [I + (1-a) * dt * D] * old_sol.

  // e -> e + 0.5 * dt * [old_source

  // Do the first solve.
  solve(solver, solver->A, solver->b, solver->x);

  // Do the second solve.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, sol2);
}

#ifdef __cplusplus
}
#endif

