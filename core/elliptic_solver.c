#include "core/elliptic_solver.h"
#include "core/hypre_helpers.h"

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
//  HYPRE_Solver pc; // Preconditioner

  bool matrix_created;

  // Index space.
  index_space_t* index_space;

  // Context pointer, virtual table.
  void* context;
  elliptic_solver_vtable vtable;
};

static void initialize(elliptic_solver_t* solver)
{
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
  solver->name = strdup(name);
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
//  if (solver->pc != NULL)
//    HYPRE_ParaSailsDestroy(solver->pc);
  HYPRE_ParCSRHybridDestroy(solver->solver);
  HYPRE_IJMatrixDestroy(solver->A);
  HYPRE_IJVectorDestroy(solver->x);
  HYPRE_IJVectorDestroy(solver->b);

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

static inline void copy_table_to_matrix(index_space_t* is, double_table_t* table, HYPRE_IJMatrix matrix)
{
}

static inline void copy_vector_to_array(index_space_t* is, HYPRE_IJVector vector, double* array)
{
  int N = is->high - is->low;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = is->low + i;
  HYPRE_IJVectorGetValues(vector, N, indices, array);
}

static inline void copy_array_to_vector(index_space_t* is, double* array, HYPRE_IJVector vector)
{
  HYPRE_IJVectorSetValuesFromArray(vector, is, array);
}

static inline void compute_op_matrix(elliptic_solver_t* solver, double_table_t* A, double t)
{
  solver->vtable.compute_operator_matrix(solver->context, A, t);
}

static inline void compute_source_vector(elliptic_solver_t* solver, double* source, double t)
{
  solver->vtable.compute_source_vector(solver->context, source, t);
}

static inline void apply_bcs(elliptic_solver_t* solver, double_table_t* A, double* b, double t)
{
  solver->vtable.apply_bcs(solver->context, A, b, t);
}

static inline void solve(elliptic_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
HYPRE_IJMatrixPrint(solver->A, "A.txt");
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

void elliptic_solver_solve(elliptic_solver_t* solver,
                           double t, double* solution)
{
  int N = solver->index_space->high - solver->index_space->low;

  // A -> operator matrix at time t.
  double_table_t* A = double_table_new();
  compute_op_matrix(solver, A, t);

  // Compute the source at time t.
  double b[N];
  compute_source_vector(solver, b, t);

  // Apply boundary conditions to the system.
  apply_bcs(solver, A, b, t);

  // Set up the linear system.
  if (!solver->matrix_created)
  {
    HYPRE_IJMatrixSetRowSizesFromTable(solver->A, solver->index_space, A);
    solver->matrix_created = true;
  }
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, A);
  copy_array_to_vector(solver->index_space, b, solver->b);

  // Solve the linear system.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to sol2.
  copy_vector_to_array(solver->index_space, solver->x, solution);

  // Clean up.
  double_table_free(A);
}

#ifdef __cplusplus
}
#endif

