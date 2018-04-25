// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BD_MATRIX_H
#define POLYMEC_BD_MATRIX_H

#include "core/polymec.h"

// This type represents a block diagonal matrix.
typedef struct bd_matrix_t bd_matrix_t;

// This returns a block diagonal matrix with a fixed block size.
bd_matrix_t* bd_matrix_new(size_t num_block_rows,
                           size_t block_size);

// This returns a block diagonal matrix with variable block sizes.
bd_matrix_t* var_bd_matrix_new(size_t num_block_rows,
                               size_t* block_sizes);

// Creates a clone of the given matrix.
bd_matrix_t* bd_matrix_clone(bd_matrix_t* matrix);

// Frees a local matrix.
void bd_matrix_free(bd_matrix_t* matrix);

// Returns the number of block rows in the matrix.
size_t bd_matrix_num_block_rows(bd_matrix_t* matrix);

// Returns the number of (individual) rows in the matrix.
size_t bd_matrix_num_rows(bd_matrix_t* matrix);

// Sets all values in the matrix to zero.
void bd_matrix_zero(bd_matrix_t* matrix);

// Sets all values in the matrix to zero.
void bd_matrix_add_identity(bd_matrix_t* matrix, real_t scale_factor);

// Returns the size of the given block row.
size_t bd_matrix_block_size(bd_matrix_t* matrix, int block_row);

// Inserts the given block into the given block row in the matrix.
void bd_matrix_insert_block(bd_matrix_t* matrix, int block_row, real_t* block);

// Adds the given block into the given block row in the matrix.
void bd_matrix_add_block(bd_matrix_t* matrix, int block_row, real_t* block);

// Simple matrix vector product matrix*vector -> product.
void bd_matrix_matvec(bd_matrix_t* matrix, real_t* vector, real_t* product);

// This returns a pointer to the storage for the given block row 
// in a block diagonal matrix.
real_t* bd_matrix_block(bd_matrix_t* matrix, int block_row);

// Prints the block diagonal matrix to the given stream.
void bd_matrix_fprintf(bd_matrix_t* matrix, FILE* stream);

// Solves a linear system A*X = B for the block diagonal matrix A and 
// vectors X, B, returning true if a solution was obtained, and false 
// otherwise.
bool solve_bd_system(bd_matrix_t* A, real_t* B, real_t* X);

#endif
