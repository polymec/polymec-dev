#ifndef POLYMEC_N_VECTOR_H
#define POLYMEC_N_VECTOR_H

#include <stdlib.h>
#include "core/polymec.h"

// An n_vector is an n-dimensional array that can be distributed 
// across several processes.
typedef struct n_vector_t n_vector_t;

// A vector destructor.
typedef void (*n_vector_dtor)(void* vec);

// Creates a new vector using the given underlying context and destructor.
n_vector_t* n_vector_new(const char* impl_name,
                         void* context,
                         n_vector_dtor dtor);

// Returns the internally-stored string with the name of the implementation
// of this vector.
char* n_vector_impl_name(n_vector_t* vector);

// Returns the context pointer for the given vector.
void* n_vector_context(n_vector_t* vector);

// Destroys the given vector.
void n_vector_free(n_vector_t* vector);

#ifdef __cplusplus
}
#endif

#endif
