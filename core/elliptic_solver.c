#include "core/elliptic_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

struct elliptic_solver_t
{
  char* name;

  // Linear system and solver.
  HYPRE_IJMatrix A;
  HYPRE_IJVector x, b;
  HYPRE_Solver solver;

  // Rows on this process.
  int ilow, ihigh;

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
  solver->vtable.create_solver(solver->context, &solver->solver);
  HYPRE_IJVectorGetLocalRange(solver->x, &solver->ilow, &solver->ihigh);
}

elliptic_solver_t* elliptic_solver_new(const char* name, 
                                       void* context,
                                       elliptic_solver_vtable vtable)
{
  ASSERT(name != NULL);
  ASSERT(vtable.create_matrix != NULL);
  ASSERT(vtable.create_vector != NULL);
  ASSERT(vtable.create_solver != NULL);
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
  HYPRE_ParCSRGMRESDestroy(&solver->solver);
  HYPRE_IJMatrixDestroy(&solver->A);
  HYPRE_IJVectorDestroy(&solver->x);
  HYPRE_IJVectorDestroy(&solver->b);

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

static inline void copy_vector_to_array(HYPRE_IJVector vector, double* array)
{
  int ilow, ihigh, N;
  HYPRE_IJVectorGetLocalRange(vector, &ilow, &ihigh);
  int N = ihigh - ilow;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = ilow + i;
  HYPRE_IJVectorGetValues(vector, N, indices, array);
}

static inline void compute_op_matrix(elliptic_solver_t* solver, HYPRE_IJMatrix A, double t)
{
  solver->vtable.compute_operator_matrix(solver->context, A, t);
}

static inline void compute_source_vector(elliptic_solver_t* solver, HYPRE_IJVector source, double t)
{
  solver->vtable.compute_source_vector(solver->context, source, t);
}

static inline void apply_bcs(elliptic_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, double t)
{
  solver->vtable.apply_bcs(solver->context, A, b, t);
}

static inline void solve(elliptic_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
  HYPRE_GMRESSetup(solver->solver, A, b, x);
  HYPRE_GMRESSolve(solver->solver, A, b, x);
}

void elliptic_solver_solve(elliptic_solver_t* solver,
                           double t, double* solution)
{
  // A -> operator matrix at time t2.
  compute_op_matrix(solver, solver->A, t);

  // Compute the source at time t2.
  compute_source_vector(solver, solver->x, t);

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
  apply_bcs(solver, solver->A, solver->b, t);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->x, solution);
}

#ifdef __cplusplus
}
#endif

