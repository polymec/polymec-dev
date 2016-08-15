// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SPARSE_LOCAL_MATRIX_H
#define POLYMEC_SPARSE_LOCAL_MATRIX_H

#include "core/local_matrix.h"
#include "core/adj_graph.h"

// This returns a sparse local matrix with a sparsity pattern given by a graph.
// The sparsity graph is consumed by the matrix.
local_matrix_t* sparse_local_matrix_new(adj_graph_t* sparsity);

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

// Initializes a set of ILU parameters with reasonable defaults:
// diagonal_pivot_threshold = 0.1
// row_perm = ILU_LARGE_DIAG_PERM
// drop_rule = ILU_DROP_BASIC | ILU_DROP_AREA
// drop_tolerance = 1e-4
// fill_factor = 10.0
// milu_variant = ILU_SILU
// fill_tolerance = 0.01
// norm = ILU_INF
ilu_params_t* ilu_params_new(void);

// This instructs the given sparse_local_matrix object to use Incomplete LU
// factorization with the given parameters.
void sparse_local_matrix_use_ilu(local_matrix_t* matrix,
                                 ilu_params_t* ilu_params);

// This instructs the given sparse_local_matrix object to use regular LU
// factorization instead of ILU.
void sparse_local_matrix_use_lu(local_matrix_t* matrix);

#endif
