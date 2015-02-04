// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCK_DIAGONAL_MATRIX_H
#define POLYMEC_BLOCK_DIAGONAL_MATRIX_H

#include "core/local_matrix.h"

// This returns a block diagonal matrix with a fixed block size.
local_matrix_t* block_diagonal_matrix_new(int num_block_rows,
                                          int block_size);

// This returns a block diagonal matrix with variable block sizes.
local_matrix_t* var_block_diagonal_matrix_new(int num_block_rows,
                                              int* block_sizes);

#endif
