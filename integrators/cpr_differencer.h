// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CPR_DIFFERENCER_H
#define POLYMEC_CPR_DIFFERENCER_H

#include "core/polymec.h"
#include "core/adj_graph.h"

// The cpr_differencer is an object that computes a Jacobian matrix for a 
// function F. The sparsity of the matrix is given by a graph, and the 
// differencer computes the matrix using the method of Curtis, Powell, and 
// Reed, in which the graph is colored so that the Jacobian can be computed 
// using finite differencing, with a minimum number of calls to F.
typedef struct cpr_differencer_t cpr_differencer_t;

// This matrix type represents a matrix that is present only on the local 
// subdomain and represents only the local portion of a Jacobian.
typedef struct local_matrix_t local_matrix_t;

// Creates a new Curtis-Powell-Reed differencer, associated with the given 
// function F *OR* F_dae and a graph, which represents the sparsity of the 
// Jacobian matrix. This constructor assumes a fixed block size.
cpr_differencer_t* cpr_differencer_new(MPI_Comm comm,
                                       int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                       int (*F_dae)(void* context, real_t, real_t* x, real_t* xdot, real_t* Fval),
                                       void* F_context,
                                       adj_graph_t* sparsity,
                                       int num_local_block_rows,
                                       int num_remote_block_rows,
                                       int block_size);

// Creates a new Curtis-Powell-Reed differencer, associated with the given 
// function F *OR* F_dae and a graph, which represents the sparsity of the 
// Jacobian matrix. This constructor allows a variable block size, 
// block_sizes[i], for the ith row of the matrix.
cpr_differencer_t* var_cpr_differencer_new(MPI_Comm comm,
                                           int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                           int (*F_dae)(void* context, real_t, real_t* x, real_t* xdot, real_t* Fval),
                                           void* F_context,
                                           adj_graph_t* sparsity,
                                           int num_local_block_rows,
                                           int num_remote_block_rows,
                                           int* block_sizes);

// Frees the differencer.
void cpr_differencer_free(cpr_differencer_t* diff);

// Uses the differencer to compute a generalized Jacobian matrix of the form
// J = alpha * I + beta * dF/dx + gamma * dF/d(xdot).
void cpr_differencer_compute(cpr_differencer_t* diff, 
                             real_t alpha, real_t beta, real_t gamma,  
                             real_t t, real_t* x, real_t* xdot,
                             local_matrix_t* matrix);

// The following interface defines a local matrix representation.
typedef struct local_matrix_vtable
{
  void (*dtor)(void* context); // Destructor
  void (*zero)(void* context); // Sets all matrix entries to zero.
  void (*add_identity)(void* context, real_t scale_factor); // Adds in scale_factor * I.
  void (*add_column_vector)(void* context,  // Adds in a scaled column vector.
                            real_t scale_factor,
                            int column,
                            real_t* column_vector);
  bool (*solve)(void* context, real_t* B, real_t* x); // Solves A*x = B.
  void (*fprintf)(void* context, FILE* stream); // Prints matrix to stream.
} local_matrix_vtable;

// This function can be used to create a new type of local matrix representation.
local_matrix_t* local_matrix_new(const char* name,
                                 void* context,
                                 local_matrix_vtable vtable);
 
// Frees a local matrix representation.
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
// solution is placed into x. Returns true if the solve is successful, 
// false otherwise.
bool local_matrix_solve(local_matrix_t* matrix, 
                        real_t* B,
                        real_t* x);

// Prints a text representation of the matrix to the given stream.
void local_matrix_fprintf(local_matrix_t* matrix, FILE* stream);

// This returns an object representing a block diagonal matrix with a 
// fixed block size.
local_matrix_t* block_diagonal_matrix_new(int num_block_rows,
                                          int block_size);

// This returns an object representing a block diagonal matrix with variable 
// block sizes.
local_matrix_t* var_block_diagonal_matrix_new(int num_block_rows,
                                              int* block_sizes);

// This returns an object representing a sparse local matrix with a
// sparsity pattern given by a graph, assuming a fixed block size.
local_matrix_t* sparse_local_matrix_new(adj_graph_t* sparsity,
                                        int num_block_rows,
                                        int block_size);

// This returns an object representing a sparse local matrix with variable 
// block sizes.
local_matrix_t* var_sparse_local_matrix_new(adj_graph_t* sparsity,
                                            int num_block_rows,
                                            int* block_sizes);

// The following types give options to control ILU preconditioners for the 
// nonlinear and time integrators.

// Specifies how Incomplete LU (ILU) solvers should permute the rows of the 
// matrix.
typedef enum
{
  ILU_NO_ROW_PERM,      // no row permutations
  ILU_LARGE_DIAG_PERM   // try to make the diagonal larger
} ilu_row_perm_t;

// Species the type of (L-) norm used to measure errors in ILU solves.
typedef enum
{
  ILU_L1,  
  ILU_L2,  
  ILU_LINF
} ilu_norm_t;

// Rules for dropping coefficients in ILU. See SuperLU documentation.
// These rules can be combined with bitwise OR.
extern const int ILU_DROP_BASIC;
extern const int ILU_DROP_PROWS;
extern const int ILU_DROP_COLUMN;
extern const int ILU_DROP_AREA;
extern const int ILU_DROP_DYNAMIC;
extern const int ILU_DROP_INTERP;

// Identifies a variant (if any) of modified ILU.
typedef enum
{
  ILU_SILU,           // no modified ILU
  ILU_MILU1,
  ILU_MILU2,
  ILU_MILU3
} ilu_variant_t;

// This container specifies a complete set of ILU parameters. Objects 
// of this type are garbage-collected.
typedef struct
{
  real_t diag_pivot_threshold;
  ilu_row_perm_t row_perm;
  int drop_rule;
  real_t drop_tolerance;
  real_t fill_factor;
  ilu_variant_t milu_variant;
  real_t fill_tolerance;
  ilu_norm_t norm;
} ilu_params_t;

// Initializes a set of ILU parameters with reasonable defaults.
ilu_params_t* ilu_params_new();

// This returns an object representing a sparse local matrix with a
// sparsity pattern given by a graph, assuming a fixed block size and
// using incomplete LU factorization for the solution of its linear system.
local_matrix_t* ilu_sparse_local_matrix_new(adj_graph_t* sparsity,
                                            int num_block_rows,
                                            int block_size,
                                            ilu_params_t* ilu_params);

// This returns an object representing a sparse local matrix with a
// sparsity pattern given by a graph, allowing a variable block size and
// using incomplete LU factorization for the solution of its linear system.
local_matrix_t* var_ilu_sparse_local_matrix_new(adj_graph_t* sparsity,
                                                int num_block_rows,
                                                int* block_sizes, 
                                                ilu_params_t* ilu_params);

#endif
