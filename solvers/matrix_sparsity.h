// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MATRIX_SPARSITY_H
#define POLYMEC_MATRIX_SPARSITY_H

#include "core/adj_graph.h"
#include "core/exchanger.h"

/// \addtogroup solvers solvers
///@{

/// \class matrix_sparsity
/// A matrix sparsity object represents the non-zero structure of a sparse 
/// matrix, including distributed matrices that are partitioned into 
/// contiguous sets of rows stored on several MPI processes. This data 
/// structure is similar in nature to an adjacency graph (adj_graph), but 
/// is defined in terms of a global index space spanning all MPI processes, 
/// instead of a separate local index space for each process.
typedef struct matrix_sparsity_t matrix_sparsity_t;

/// Creates a new matrix sparsity pattern on the given MPI communicator
/// with the given row_distribution array, which contains nprocs+1 elements. 
/// Specifically, row_distribution[p] gives the global index of the first 
/// row stored on process p, and the number of rows local to process p can 
/// be calculated as row_distribution[p+1] - row_distribution[p].
/// \memberof matrix_sparsity
matrix_sparsity_t* matrix_sparsity_new(MPI_Comm comm, 
                                       index_t* row_distribution);

/// Creates a new matrix sparsity pattern using the given sparsity pattern, 
/// and interpreting each non-zero element as a block of the given size. The 
/// resulting sparsity pattern will contain non-zero elements for each of those
/// within the blocks of the originally given sparsity pattern.
/// \memberof matrix_sparsity
matrix_sparsity_t* matrix_sparsity_with_block_size(matrix_sparsity_t* sparsity,
                                                   size_t block_size);

/// Creates a new matrix sparsity pattern using the given sparsity pattern, 
/// and interpreting each non-zero element as a block of size specified for 
/// each row. The resulting sparsity pattern will contain non-zero elements for 
/// each of those within the blocks of the originally given sparsity pattern.
/// Here, block_sizes[i] contains the block size of the ith row in the original
/// sparsity pattern.
/// \memberof matrix_sparsity
matrix_sparsity_t* matrix_sparsity_with_block_sizes(matrix_sparsity_t* sparsity,
                                                    size_t* block_sizes);

/// Creates a new sparsity pattern using the given adjacency graph (which 
/// provides a topology connecting local indices) and an exchanger (which 
/// can be used to construct a mapping from local to global indices for 
/// graph edges that cross MPI process boundaries). The sparsity pattern 
/// is defined on the same MPI communicator as the adjacency graph.
/// \note This operation involves communication if the number of processes
/// in the graph's communicator is greater than 1. On single processes,
/// it's okay to pass NULL to the exchanger.
/// \memberof matrix_sparsity
matrix_sparsity_t* matrix_sparsity_from_graph(adj_graph_t* graph,
                                              exchanger_t* ex);

/// Creates a clone (deep copy) of the given matrix sparsity pattern.
/// \memberof matrix_sparsity
matrix_sparsity_t* matrix_sparsity_clone(matrix_sparsity_t* sparsity);

/// Frees the given matrix sparsity pattern.
/// \memberof matrix_sparsity
void matrix_sparsity_free(matrix_sparsity_t* sparsity);

/// Returns the communicator for this sparsity object.
/// \memberof matrix_sparsity
MPI_Comm matrix_sparsity_comm(matrix_sparsity_t* sparsity);

/// Returns the number of global rows in the sparsity object.
/// \memberof matrix_sparsity
size_t matrix_sparsity_num_global_rows(matrix_sparsity_t* sparsity);

/// Returns the number of local rows in the sparsity object.
/// \memberof matrix_sparsity
size_t matrix_sparsity_num_local_rows(matrix_sparsity_t* sparsity);

/// Returns an internal pointer to an array that describes the distribution
/// of rows on each process. The pth entry in the array holds the global 
/// index of the first row that is stored on process p, and the number of 
/// rows locally stored on process p is the difference of the (p+1)th and 
/// pth entries.
/// \memberof matrix_sparsity
index_t* matrix_sparsity_row_distribution(matrix_sparsity_t* sparsity);

/// Returns the total number of nonzero entries in the sparsity pattern.
/// \memberof matrix_sparsity
size_t matrix_sparsity_num_nonzeros(matrix_sparsity_t* sparsity);

/// Sets the number of columns for the given row in the sparsity object.
/// \memberof matrix_sparsity
void matrix_sparsity_set_num_columns(matrix_sparsity_t* sparsity, 
                                     index_t row, 
                                     size_t num_columns);

/// Returns the number of columns for the given row in the sparsity object.
/// \memberof matrix_sparsity
size_t matrix_sparsity_num_columns(matrix_sparsity_t* sparsity, 
                                   index_t row);

/// Returns an internal pointer to the array of column indices for the given
/// row of the sparsity object. This can be used to retrieve or to set the 
/// columns for the row, but matrix_sparsity_set_num_columns must be called 
/// before setting the column indices for a row.
/// \memberof matrix_sparsity
index_t* matrix_sparsity_columns(matrix_sparsity_t* sparsity, index_t row);

/// Returns true if the given row in the sparsity pattern includes the given 
/// column, false if not.
/// \memberof matrix_sparsity
bool matrix_sparsity_contains(matrix_sparsity_t* sparsity, 
                              index_t row, 
                              index_t column);

/// Allows iteration over the (global) row indices for the rows stored on the 
/// local MPI process in the sparsity. Returns true if more rows are found, 
/// false if not. Set pos to 0 to reset an iteration.
/// \memberof matrix_sparsity
bool matrix_sparsity_next_row(matrix_sparsity_t* sparsity, 
                              int* pos, 
                              index_t* row);

/// Allows iteration over the column indices for the given row in the sparsity.
/// Returns true if more columns are found, false if not. Set pos to 0 to 
/// reset an iteration.
/// \memberof matrix_sparsity
bool matrix_sparsity_next_column(matrix_sparsity_t* sparsity, 
                                 index_t row,
                                 int* pos, 
                                 index_t* column);

/// Prints a textual representation of the matrix sparsity to the given file.
/// \memberof matrix_sparsity
void matrix_sparsity_fprintf(matrix_sparsity_t* sparsity, FILE* stream);

///@}

#endif
