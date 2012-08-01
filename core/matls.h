#ifndef ARBI_MATLS_H
#define ARBI_MATLS_H

#include "core/arbi.h"
#include "core/mat.h"
#include "core/vec.h"

#ifdef __cplusplus
extern "C" {
#endif

// matls is a matrix-based sparse linear solver class.
typedef struct matls_t matls_t;

// A function pointer type for constructing a vector.
typedef vec_t* (*matls_vector_func)(void*, int);

// A function pointer type for constructing a matrix.
typedef mat_t* (*matls_matrix_func)(void*, int);

// A function pointer type for solving a linear system.
typedef int (*matls_solve_func)(void*, mat_t*, vec_t*, vec_t*);

// A destructor function for the context object (if any).
typedef void (*matls_dtor)(void*);

// This virtual table must be implemented by any matls subclass.
typedef struct 
{
  matls_vector_func             vector;
  matls_matrix_func             matrix;
  matls_solve_func              solve;
  matls_dtor                    dtor;
} matls_vtable;

// Constructs a matrix-based linear solver.
matls_t* matls_new(void* context, const char* name, matls_vtable vtable);

// Creates and returns an Nx1 vector suitable for use by this linear solver.
vec_t* matls_vector(matls_t* solver, int N);

// Creates and returns an NxN matrix suitable for use by this linear solver.
mat_t* matls_matrix(matls_t* solver, int N);

// Destroys the given (matrix-based) linear solver.
void matls_free(matls_t* solver);

// Solves the linear system Ax = b using the given linear solver.
void matls_solve(matls_t* solver, mat_t* A, vec_t* x, vec_t* b);

#ifdef __cplusplus
}
#endif

#endif

