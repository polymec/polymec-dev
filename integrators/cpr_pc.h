// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CPR_PC_H
#define POLYMEC_CPR_PC_H

#include "core/sparse_local_matrix.h"
#include "integrators/newton_pc.h"

// The Curtis-Powell-Reed preconditioner is a Newton preconditioner that uses
// the method of Curtis, Powell and Reed to approximate the entries of a 
// Jacobian matrix automatically, given only the residual (or right-hand side) 
// function.

// Each of the following creates a Curtis-Powell-Reed preconditioner for the 
// given function F(t, x), to construct matrices with a sparsity pattern 
// expressed by the given sparsity graph, number of block rows, and block size. 
// The nature of the sparsity graph (i.e. whether it is a block graph or not) is 
// inferred from the number of block rows and the block size. 

// Block Jacobi, fixed block size.
newton_pc_t* block_jacobi_cpr_pc_from_function(MPI_Comm comm,
                                               void* context,
                                               int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                               void (*dtor)(void* context),
                                               adj_graph_t* sparsity,
                                               int num_local_block_rows,
                                               int num_remote_block_rows,
                                               int block_size);
                                        
// Block Jacobi, variable block size.
newton_pc_t* var_block_jacobi_cpr_pc_from_function(MPI_Comm comm,
                                                   void* context,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                   void (*dtor)(void* context),
                                                   adj_graph_t* sparsity,
                                                   int num_local_block_rows,
                                                   int num_remote_block_rows,
                                                   int* block_sizes);
 
// LU decomposition, fixed block size.
newton_pc_t* lu_cpr_pc_from_function(MPI_Comm comm,
                                     void* context,
                                     int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                     void (*dtor)(void* context),
                                     adj_graph_t* sparsity,
                                     int num_local_block_rows,
                                     int num_remote_block_rows,
                                     int block_size);
                                        
// LU decomposition, variable block size.
newton_pc_t* var_lu_cpr_pc_from_function(MPI_Comm comm,
                                         void* context,
                                         int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                         void (*dtor)(void* context),
                                         adj_graph_t* sparsity,
                                         int num_local_block_rows,
                                         int num_remote_block_rows,
                                         int* block_sizes);
 
// ILU decomposition, fixed block size.
newton_pc_t* ilu_cpr_pc_from_function(MPI_Comm comm,
                                      void* context,
                                      int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                      void (*dtor)(void* context),
                                      adj_graph_t* sparsity,
                                      int num_local_block_rows,
                                      int num_remote_block_rows,
                                      int block_size,
                                      ilu_params_t* ilu_params);
                                        
// ILU decomposition, variable block size.
newton_pc_t* var_ilu_cpr_pc_from_function(MPI_Comm comm,
                                          void* context,
                                          int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                          void (*dtor)(void* context),
                                          adj_graph_t* sparsity,
                                          int num_local_block_rows,
                                          int num_remote_block_rows,
                                          int* block_sizes,
                                          ilu_params_t* ilu_params);
 
// Each of the following creates a Curtis-Powell-Reed preconditioner for the 
// given differential algebraic equation (DAE) function F(t, x, xdot), to 
// construct matrices with a sparsity pattern expressed by the given sparsity 
// graph, number of block rows, and block size. The nature of the sparsity 
// graph (i.e. whether it is a block graph or not) is inferred from the number 
// of block rows and the block size. 

// Block Jacobi, fixed block size.
newton_pc_t* block_jacobi_cpr_pc_from_dae_function(MPI_Comm comm,
                                                   void* context,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                   void (*dtor)(void* context),
                                                   adj_graph_t* sparsity,
                                                   int num_local_block_rows,
                                                   int num_remote_block_rows,
                                                   int block_size);

// Block Jacobi, variable block size.
newton_pc_t* var_block_jacobi_cpr_pc_from_dae_function(MPI_Comm comm,
                                                       void* context,
                                                       int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                       void (*dtor)(void* context),
                                                       adj_graph_t* sparsity,
                                                       int num_local_block_rows,
                                                       int num_remote_block_rows,
                                                       int* block_sizes);

// LU decomposition, fixed block size.
newton_pc_t* lu_cpr_pc_from_dae_function(MPI_Comm comm,
                                         void* context,
                                         int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                         void (*dtor)(void* context),
                                         adj_graph_t* sparsity,
                                         int num_local_block_rows,
                                         int num_remote_block_rows,
                                         int block_size);

// LU decomposition, variable block size.
newton_pc_t* var_lu_cpr_pc_from_dae_function(MPI_Comm comm,
                                             void* context,
                                             int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                             void (*dtor)(void* context),
                                             adj_graph_t* sparsity,
                                             int num_local_block_rows,
                                             int num_remote_block_rows,
                                             int* block_sizes);

// ILU decomposition, fixed block size.
newton_pc_t* ilu_cpr_pc_from_dae_function(MPI_Comm comm,
                                          void* context,
                                          int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                          void (*dtor)(void* context),
                                          adj_graph_t* sparsity,
                                          int num_local_block_rows,
                                          int num_remote_block_rows,
                                          int block_size,
                                          ilu_params_t* ilu_params);

// ILU decomposition, variable block size.
newton_pc_t* var_ilu_cpr_pc_from_dae_function(MPI_Comm comm,
                                              void* context,
                                              int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                              void (*dtor)(void* context),
                                              adj_graph_t* sparsity,
                                              int num_local_block_rows,
                                              int num_remote_block_rows,
                                              int* block_sizes,
                                              ilu_params_t* ilu_params);

#endif

