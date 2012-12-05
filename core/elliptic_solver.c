#include "core/elliptic_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

struct elliptic_solver_t
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
  elliptic_solver_vtable vtable;
};

static void initialize(elliptic_solver_t* solver)
{
  solver->vtable.create_matrix(solver->context, &solver->A);
  solver->vtable.create_vector(solver->context, &solver->x);
  solver->vtable.create_vector(solver->context, &solver->b);
  solver->vtable.create_ksp(solver->context, &solver->ksp);
}

elliptic_solver_t* elliptic_solver_new(const char* name, 
                                         void* context,
                                         elliptic_solver_vtable vtable)
{
  ASSERT(name != NULL);
  ASSERT(vtable.create_matrix != NULL);
  ASSERT(vtable.create_vector != NULL);
  ASSERT(vtable.create_ksp != NULL);
  ASSERT(vtable.compute_operator_matrix != NULL);
  ASSERT(vtable.apply_bcs != NULL);

  elliptic_solver_t* solver = malloc(sizeof(elliptic_solver_t));
  solver->name = strdup(name);
  solver->context = context;
  solver->vtable = vtable;

  // Make sure the solver is initialized.
  initialize(solver);

  return solver;
}

void elliptic_solver_free(elliptic_solver_t* solver)
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

char* elliptic_solver_name(elliptic_solver_t* solver)
{
  return solver->name;
}

void* elliptic_solver_context(elliptic_solver_t* solver)
{
  return solver->context;
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

static inline void compute_op_matrix(elliptic_solver_t* solver, Mat A, double t)
{
  solver->vtable.compute_operator_matrix(solver->context, A, t);
}

static inline void compute_source_vector(elliptic_solver_t* solver, Vec source, double t)
{
  solver->vtable.compute_source_vector(solver->context, source, t);
}

static inline void apply_bcs(elliptic_solver_t* solver, Mat A, Vec b, double t)
{
  solver->vtable.apply_bcs(solver->context, A, b, t);
}

static inline void solve(elliptic_solver_t* solver, Mat A, Vec b, Vec x)
{
  KSPSetOperators(solver->ksp, A, A, SAME_NONZERO_PATTERN);
  KSPSolve(solver->ksp, b, x);
}

void elliptic_solver_solve(elliptic_solver_t* solver,
                           double t, double* solution)
{
  // A -> operator matrix at time t2.
  compute_op_matrix(solver, solver->A, t);

  // Compute the source at time t2.
  compute_source_vector(solver, solver->x, t);

  // Apply boundary conditions to the system.
  VecSet(solver->b, 0.0);
  apply_bcs(solver, solver->A, solver->b, t);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, solution);
}

#ifdef __cplusplus
}
#endif

