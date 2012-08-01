#include <stdlib.h>
#include <string.h>
#include "core/matls.h"

#ifdef __cplusplus
extern "C" {
#endif

// matls is a matrix-based sparse linear solver class.
struct matls_t 
{
  char*        name;
  void*        context;
  matls_vtable vtable;
};

matls_t* matls_new(void* context, const char* name, matls_vtable vtable)
{
  ASSERT(vtable.vector != NULL);
  ASSERT(vtable.matrix != NULL);
  ASSERT(vtable.solve != NULL);
  matls_t* ls = malloc(sizeof(matls_t));
  ls->context = context;
  ls->name = strdup(name);
  ls->vtable = vtable;
  return ls;
}

vec_t* matls_vector(matls_t* solver, int N)
{
  return solver->vtable.vector(solver->context, N);
}

mat_t* matls_matrix(matls_t* solver, int N)
{
  return solver->vtable.matrix(solver->context, N);
}

void matls_free(matls_t* solver)
{
  free(solver->name);
  if (solver->vtable.dtor != NULL)
    solver->vtable.dtor(solver->context);
}

void matls_solve(matls_t* solver, mat_t* A, vec_t* x, vec_t* b)
{
  solver->vtable.solve(solver->context, A, x, b);
}

#ifdef __cplusplus
}
#endif

