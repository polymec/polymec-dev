// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LOCAL_MATRIX_H
#define POLYMEC_LOCAL_MATRIX_H

#include "core/polymec.h"

// This matrix type represents a matrix that is present only on the local 
// subdomain, plus a built-in solver. It is intended only for solving 
// EXCEEDINGLY SIMPLE OR SMALL linear systems that may have sparsity or 
// other characteristics that would make simple LAPACK solves inconvenient.
typedef struct local_matrix_t local_matrix_t;

// The following interface defines a local matrix representation.
// Use this with local_matrix_new to construct a new type of local matrix.
typedef struct local_matrix_vtable
{
  void (*dtor)(void* context); // Destructor
  void (*zero)(void* context); // Sets all matrix entries to zero.
  void (*add_identity)(void* context, real_t scale_factor); // Adds in scale_factor * I.
  void (*add_column_vector)(void* context,  // Adds in a scaled column vector.
                            real_t scale_factor,
                            int column,
                            real_t* column_vector);
  bool (*solve)(void* context, real_t* B); // Solves A*x = B in place.
  void (*fprintf)(void* context, FILE* stream); // Prints matrix to stream.
  real_t (*value)(void* context, int i, int j); // A(i, j)
  void (*set_value)(void* context, int i, int j, real_t value); // A(i, j) = value
} local_matrix_vtable;

// This can be used to create a new type of local matrix representation.
local_matrix_t* local_matrix_new(const char* name,
                                 void* context,
                                 local_matrix_vtable vtable);
 
// Frees a local matrix.
void local_matrix_free(local_matrix_t* matrix);

// Returns an internally stored string containing the name of the local 
// matrix representation.
char* local_matrix_name(local_matrix_t* matrix);

// Sets all entries in a local matrix to zero.
void local_matrix_zero(local_matrix_t* matrix);

// Adds a scaled column vector into a local matrix.
void local_matrix_add_column_vector(local_matrix_t* matrix, 
                                    real_t scale_factor,
                                    int column,
                                    real_t* column_vector);

// Adds a scaled identity matrix into a local matrix.
void local_matrix_add_identity(local_matrix_t* matrix, 
                               real_t scale_factor);

// Solves the linear system A*x = B, where A is the local matrix. The 
// solution is placed into B. Returns true if the solve is successful, 
// false otherwise.
bool local_matrix_solve(local_matrix_t* matrix, 
                        real_t* B);

// Prints a text representation of the matrix to the given stream.
void local_matrix_fprintf(local_matrix_t* matrix, FILE* stream);

// Returns the value of the matrix at the ith row, jth column. This 
// is useful for debugging, but probably not for doing heavy lifting.
real_t local_matrix_value(local_matrix_t* matrix, int i, int j);

// Sets the value of the matrix at the ith row, jth column. This 
// is useful for fiddling, but probably not for doing heavy lifting.
// If (i, j) is not a non-zero matrix index, this call has no effect.
void local_matrix_set_value(local_matrix_t* matrix, int i, int j, real_t value);

#endif
