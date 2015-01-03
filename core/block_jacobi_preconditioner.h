// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCK_JACOBI_PRECONDITIONER_H
#define POLYMEC_BLOCK_JACOBI_PRECONDITIONER_H

#include "core/preconditioner.h"
#include "core/adj_graph.h"

// Creates a block Jacobi preconditioner from a rule that computes the 
// (block) diagonal. The diagonal D is stored in block-major order. Specifically, 
// it is an array consisting of the components of column-major dense matrices 
// stacked end-to-end.
preconditioner_t* block_jacobi_preconditioner_new(void* context,
                                                  void (*compute_diagonal)(void* context, int block_size, real_t* D),
                                                  void (*dtor)(void* context),
                                                  int num_block_rows,
                                                  int block_size);

// Creates a block Jacobi preconditioner that uses data in the given 
// (block-major) array as its diagonal. The preconditioner does NOT manage the 
// array -- this array is assumed to be valid for the lifetime of the preconditioner.
preconditioner_t* block_jacobi_preconditioner_from_array(real_t* array,
                                                         int num_block_rows,
                                                         int block_size);

#endif

