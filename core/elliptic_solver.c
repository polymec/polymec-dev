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

  // Index space.
  index_space_t* index_space;

  // Context pointer, virtual table.
  void* context;
  elliptic_solver_vtable vtable;
};

static void initialize(elliptic_solver_t* solver)
{
  MPI_Comm comm = solver->index_space->comm;
  int low = solver->index_space->low;
  int high = solver->index_space->high;
  HYPRE_IJMatrixCreate(comm, low, high, low, high, &solver->A);
  HYPRE_IJMatrixSetObjectType(solver->A, HYPRE_PARCSR);
  HYPRE_IJVectorCreate(comm, low, high, &solver->x);
  HYPRE_IJVectorSetObjectType(solver->x, HYPRE_PARCSR);
  HYPRE_IJVectorCreate(comm, low, high, &solver->b);
  HYPRE_IJVectorSetObjectType(solver->b, HYPRE_PARCSR);
  HYPRE_ParCSRGMRESCreate(comm, &solver->solver);
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

  // Make sure the solver is initialized.
  initialize(solver);

  return solver;
}

void elliptic_solver_free(elliptic_solver_t* solver)
{
  HYPRE_ParCSRGMRESDestroy(solver->solver);
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
  // FIXME: Pre-allocate matrix entries if it hasn't been done already.
#if 0
  HYPRE_IJMatrixSetRowSizesFromTable(matrix, table);
#endif
  HYPRE_IJMatrixSetValuesFromTable(matrix, is, table);
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
  HYPRE_ParCSRMatrix Aobj;
  HYPRE_IJMatrixGetObject(solver->A, (void**)&Aobj);
  HYPRE_ParVector xobj, bobj;
  HYPRE_IJVectorGetObject(solver->x, (void**)&xobj);
  HYPRE_IJVectorGetObject(solver->b, (void**)&bobj);

  HYPRE_GMRESSetup(solver->solver, (HYPRE_Matrix)Aobj, (HYPRE_Vector)bobj, (HYPRE_Vector)xobj);
  HYPRE_GMRESSolve(solver->solver, (HYPRE_Matrix)Aobj, (HYPRE_Vector)bobj, (HYPRE_Vector)xobj);
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
  copy_table_to_matrix(solver->index_space, A, solver->A);
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

