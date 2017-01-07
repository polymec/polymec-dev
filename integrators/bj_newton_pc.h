// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCK_JACOBI_NEWTON_PC_H
#define POLYMEC_BLOCK_JACOBI_NEWTON_PC_H

#include "core/adj_graph.h"
#include "integrators/newton_pc.h"

// The block Jacobi Newton PC implements block Jacobi preconditioning given 
// a function that computes the diagonal block of the Jacobian.

// Creates a block Jacobi Newton preconditioner using the given function 
// compute_diag_block, which computes the diagonal block of the Jacobian 
// matrix alpha*I + beta*dFdx + gamma*dFd(x_dot) for the ith block row at 
// time t and given x and x_dot. The diagonal block is stored in column-major 
// format (so it is compatible with calls to LAPACK).
newton_pc_t* bj_newton_pc_new(void* context,
                              void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                              void (*dtor)(void* context),
                              newton_pc_side_t side,
                              int num_block_rows,
                              int block_size);
 
// Creates a block Jacobi Newton preconditioner with variable block sizes.
newton_pc_t* var_bj_newton_pc_new(void* context,
                                  void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                                  void (*dtor)(void* context),
                                  newton_pc_side_t side,
                                  int num_block_rows,
                                  int* block_sizes);

// Creates a block Jacobi Newton preconditioner that uses Curtis, Powell, and 
// Reed's method (On the estimation of sparse Jacobian Matrices, 
// J. Inst. Math. Appl., 13 (1974), pp. 117-119) for computing the Jacobian 
// matrix alpha*I + beta*dFdx using finite differences with graph coloring to 
// minimize function calls.
newton_pc_t* cpr_bj_newton_pc_new(MPI_Comm comm,
                                  void* context,
                                  int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                  void (*dtor)(void* context),
                                  newton_pc_side_t side,
                                  adj_graph_t* sparsity,
                                  int num_local_block_rows,
                                  int num_remote_block_rows,
                                  int block_size);

// This is the variable-block-size version of the Curtis-Powell-Reed block 
// Jacobi preconditioner described above.
newton_pc_t* var_cpr_bj_newton_pc_new(MPI_Comm comm,
                                      void* context,
                                      int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                      void (*dtor)(void* context),
                                      newton_pc_side_t side,
                                      adj_graph_t* sparsity,
                                      int num_local_block_rows,
                                      int num_remote_block_rows,
                                      int* block_sizes);

// Creates a block Jacobi Newton preconditioner that uses Curtis, Powell, and 
// Reed's method (On the estimation of sparse Jacobian Matrices, 
// J. Inst. Math. Appl., 13 (1974), pp. 117-119) for computing the Jacobian 
// matrix alpha*I + beta*dFdx + gamma*dFdxdot appropriate for Differential 
// Algebraic Equaions, using finite differences with graph coloring to 
// minimize function calls.
newton_pc_t* dae_cpr_bj_newton_pc_new(MPI_Comm comm,
                                      void* context,
                                      int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                      void (*dtor)(void* context),
                                      adj_graph_t* sparsity,
                                      int num_local_block_rows,
                                      int num_remote_block_rows,
                                      int block_size);

// This is the variable-block-size version of the Differential Algebraic 
// Equation (DAE) Curtis-Powell-Reed block Jacobi preconditioner described above.
newton_pc_t* var_dae_cpr_bj_newton_pc_new(MPI_Comm comm,
                                          void* context,
                                          int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                          void (*dtor)(void* context),
                                          adj_graph_t* sparsity,
                                          int num_local_block_rows,
                                          int num_remote_block_rows,
                                          int* block_sizes);

#endif

