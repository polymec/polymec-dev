#include "core/polymec.h"
#include "core/sparse_lin_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

struct sparse_lin_solver_t
{
  char* name;
  void* context;
  sparse_lin_solver_vtable vtable;

  // Solver metadata.
  double res_norm;
  int nli, nps;
};

sparse_lin_solver_t* sparse_lin_solver_new(const char* name, void* context,
                                           sparse_lin_solver_vtable vtable)
{
  ASSERT(vtable.solve != NULL);
  sparse_lin_solver_t* solver = malloc(sizeof(sparse_lin_solver_t));
  solver->name = strdup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->res_norm = 0.0;
  solver->nli = 0;
  solver->nps = 0;
  return solver;
}

void sparse_lin_solver_free(sparse_lin_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  free(solver->name);
  free(solver);
}

void sparse_lin_solver_solve(sparse_lin_solver_t* solver, N_Vector x, N_Vector b)
{
  solver->res_norm = 0.0;
  solver->nli = 0;
  solver->nps = 0;
  solver->vtable.solve(solver->context, x, b, &solver->res_norm, &solver->nli, &solver->nps);
}

void sparse_lin_solver_get_info(sparse_lin_solver_t* solver,
                                double* res_l2_norm,
                                int* num_linear_iterations,
                                int* num_precond_solves)
{
  *res_l2_norm = solver->res_norm;
  *num_linear_iterations = solver->nli;
  *num_precond_solves = solver->nps;
}

#ifdef __cplusplus
}
#endif

