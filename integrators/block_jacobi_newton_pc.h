// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCK_JACOBI_NEWTON_PC_H
#define POLYMEC_BLOCK_JACOBI_NEWTON_PC_H

#include "integrators/newton_pc.h"

// The block Jacobi Newton PC implements block Jacobi preconditioning given 
// a function that computes the diagonal block of the Jacobian.

// Creates a block Jacobi Newton preconditioner using the given function 
// compute_diag_block, which computes the diagonal block of the Jacobian 
// matrix alpha*I + beta*dFdx + gamma*dFd(x_dot) for the ith block row at 
// time t and given x and x_dot. The diagonal block is stored in column-major 
// format (so it is compatible with calls to LAPACK).
newton_pc_t* block_jacobi_newton_pc_new(void* context,
                                        void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                                        void (*dtor)(void* context),
                                        newton_pc_side_t side,
                                        int num_block_rows,
                                        int block_size);
 
 // Creates a block Jacobi Newton preconditioner with variable block sizes.
newton_pc_t* var_block_jacobi_newton_pc_new(void* context,
                                            void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                                            void (*dtor)(void* context),
                                            newton_pc_side_t side,
                                            int num_block_rows,
                                            int* block_sizes);

#endif

