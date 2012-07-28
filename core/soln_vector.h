#ifndef ARBI_SOLN_VECTOR_H
#define ARBI_SOLN_VECTOR_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// The representation of a numerical solution vector.
typedef struct 
{
  // Total number of elements in the solution vector.
  int size;
  // Number of components.
  int components;
  // A pointer to the linearized array data.
  real_t* data;
} soln_vector_t;

// Construct a new vector containing zeros.
soln_vector_t* soln_vector_new(int size, int components);

// Create a copy of the given vector.
soln_vector_t* soln_vector_clone(soln_vector_t* soln);

// Destroy the solution vector.
void soln_vector_free(soln_vector_t* soln);

#ifdef __cplusplus
}
#endif

#endif

