#ifndef POLYMEC_N_VECTOR_H
#define POLYMEC_N_VECTOR_H

#include <stdlib.h>
#include "core/polymec.h"

// An n_vector is an n-dimensional array that can be distributed 
// across several processes.
typedef struct n_vector_t n_vector_t;

// A vector clone operation.
typedef void* (*n_vector_clone_func)(void* vec);

// A data access function.
typedef double* (*n_vector_data_func)(void* vec);

// A function to return the vector's local size.
typedef int (*n_vector_local_size_func)(void* vec);

// A function to return the vector's global size.
typedef int (*n_vector_global_size_func)(void* vec);

// A vector destructor.
typedef void (*n_vector_dtor)(void* vec);

// Any implementation of a vector must provide this virtual table.
typedef struct
{
  n_vector_clone_func clone;
  n_vector_data_func data;
  n_vector_local_size_func local_size;
  n_vector_global_size_func global_size;
  n_vector_dtor dtor;
} n_vector_vtable;

// Creates a new vector using the given underlying implementation and 
// context.
n_vector_t* n_vector_new(const char* impl_name,
                         void* context,
                         n_vector_vtable vtable);

// Returns the internally-stored string with the name of the implementation
// of this vector.
char* n_vector_impl_name(n_vector_t* vector);

// Returns the context pointer for the given vector.
void* n_vector_context(n_vector_t* vector);

// Clones the given vector, returning a fresh copy.
n_vector_t* n_vector_clone(n_vector_t* vector);

// Returns a pointer to the locally-stored array of data within the 
// given vector.
double* n_vector_data(n_vector_t* vector);

// Returns the size of the locally-stored data in the vector.
int n_vector_local_size(n_vector_t* vector);

// Returns the global size of the data in the vector.
int n_vector_global_size(n_vector_t* vector);

// Destroys the given vector.
void n_vector_free(n_vector_t* vector);

#ifdef __cplusplus
}
#endif

#endif
