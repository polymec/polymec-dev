// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/block_diagonal_matrix.h"
#include "core/sparse_local_matrix.h"
#include "integrators/cpr_pc.h"
#include "integrators/cpr_differencer.h"

// This preconditioner consists of a differencer and a given matrix.
typedef struct 
{
  cpr_differencer_t* diff;
  local_matrix_t* P;
} cpr_pc_t;

static void cpr_pc_compute_p(void* context, 
                             real_t alpha, real_t beta, real_t gamma, 
                             real_t t, real_t* x, real_t* xdot)
{
  cpr_pc_t* pc = context;
  cpr_differencer_compute(pc->diff, alpha, beta, gamma, t, x, xdot, pc->P);
}

static bool cpr_pc_solve(void* context, real_t* R)
{
  cpr_pc_t* pc = context;
  return local_matrix_solve(pc->P, R);
}

static void cpr_pc_dtor(void* context)
{
  cpr_pc_t* pc = context;
  cpr_differencer_free(pc->diff);
  local_matrix_free(pc->P);
  polymec_free(pc);
}

static newton_pc_t* cpr_pc_from_function_and_matrix(MPI_Comm comm,
                                                    int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                    int (*dae_F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                    void* context,
                                                    void (*dtor)(void* context),
                                                    adj_graph_t* sparsity,
                                                    int num_local_block_rows,
                                                    int num_remote_block_rows,
                                                    int block_size,
                                                    local_matrix_t* P)
{
  cpr_pc_t* pc = polymec_malloc(sizeof(cpr_pc_t));
  pc->diff = cpr_differencer_new(comm, F, NULL, context, dtor,
                                 sparsity, num_local_block_rows,
                                 num_remote_block_rows, block_size);
  pc->P = P;
  newton_pc_vtable vtable = {.compute_p = cpr_pc_compute_p,
                             .solve = cpr_pc_solve,
                             .dtor = cpr_pc_dtor};
  return newton_pc_new("Curtis-Powell-Reed preconditioner", pc, vtable);
}

static newton_pc_t* var_cpr_pc_from_function_and_matrix(MPI_Comm comm,
                                                        int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                        int (*dae_F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                        void* context,
                                                        void (*dtor)(void* context),
                                                        adj_graph_t* sparsity,
                                                        int num_local_block_rows,
                                                        int num_remote_block_rows,
                                                        int* block_sizes,
                                                        local_matrix_t* P)
{
  cpr_pc_t* pc = polymec_malloc(sizeof(cpr_pc_t));
  pc->diff = var_cpr_differencer_new(comm, F, NULL, context, dtor,
                                     sparsity, num_local_block_rows,
                                     num_remote_block_rows, block_sizes);
  pc->P = P;
  newton_pc_vtable vtable = {.compute_p = cpr_pc_compute_p,
                             .solve = cpr_pc_solve,
                             .dtor = cpr_pc_dtor};
  return newton_pc_new("Curtis-Powell-Reed preconditioner", pc, vtable);
}

newton_pc_t* block_jacobi_cpr_pc_from_function(MPI_Comm comm,
                                               int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                               void* context,
                                               void (*dtor)(void* context),
                                               adj_graph_t* sparsity,
                                               int num_local_block_rows,
                                               int num_remote_block_rows,
                                               int block_size)
{
  local_matrix_t* P = block_diagonal_matrix_new(num_local_block_rows, block_size);
  return cpr_pc_from_function_and_matrix(comm, F, NULL, context, dtor, 
                                         sparsity, num_local_block_rows, 
                                         num_remote_block_rows, block_size, P);
}
                                        
newton_pc_t* var_block_jacobi_cpr_pc_from_function(MPI_Comm comm,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                   void* context,
                                                   void (*dtor)(void* context),
                                                   adj_graph_t* sparsity,
                                                   int num_local_block_rows,
                                                   int num_remote_block_rows,
                                                   int* block_sizes)
{
  local_matrix_t* P = var_block_diagonal_matrix_new(num_local_block_rows, block_sizes);
  return var_cpr_pc_from_function_and_matrix(comm, F, NULL, context, dtor, 
                                             sparsity, num_local_block_rows, 
                                             num_remote_block_rows, block_sizes, P);
}
 
newton_pc_t* lu_cpr_pc_from_function(MPI_Comm comm,
                                     int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                     void* context,
                                     void (*dtor)(void* context),
                                     adj_graph_t* sparsity,
                                     int num_local_block_rows,
                                     int num_remote_block_rows,
                                     int block_size)
{
  adj_graph_t* block_sparsity = adj_graph_new_with_block_size(sparsity, block_size);
  local_matrix_t* P = sparse_local_matrix_new(block_sparsity);
  adj_graph_free(block_sparsity);
  return cpr_pc_from_function_and_matrix(comm, F, NULL, context, dtor, 
                                         sparsity, num_local_block_rows, 
                                         num_remote_block_rows, block_size, P);
}
                                        
newton_pc_t* var_lu_cpr_pc_from_function(MPI_Comm comm,
                                         int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                         void* context,
                                         void (*dtor)(void* context),
                                         adj_graph_t* sparsity,
                                         int num_local_block_rows,
                                         int num_remote_block_rows,
                                         int* block_sizes)
{
  adj_graph_t* block_sparsity = adj_graph_new_with_block_sizes(sparsity, block_sizes);
  local_matrix_t* P = sparse_local_matrix_new(block_sparsity);
  adj_graph_free(block_sparsity);
  return var_cpr_pc_from_function_and_matrix(comm, F, NULL, context, dtor, 
                                             sparsity, num_local_block_rows, 
                                             num_remote_block_rows, block_sizes, P);
}
 
newton_pc_t* ilu_cpr_pc_from_function(MPI_Comm comm,
                                      int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                      void* context,
                                      void (*dtor)(void* context),
                                      adj_graph_t* sparsity,
                                      int num_local_block_rows,
                                      int num_remote_block_rows,
                                      int block_size,
                                      ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);
  adj_graph_t* block_sparsity = adj_graph_new_with_block_size(sparsity, block_size);
  local_matrix_t* P = ilu_sparse_local_matrix_new(block_sparsity, ilu_params);
  adj_graph_free(block_sparsity);
  return cpr_pc_from_function_and_matrix(comm, F, NULL, context, dtor, 
                                         sparsity, num_local_block_rows, 
                                         num_remote_block_rows, block_size, P);
}
                                        
newton_pc_t* var_ilu_cpr_pc_from_function(MPI_Comm comm,
                                          int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                          void* context,
                                          void (*dtor)(void* context),
                                          adj_graph_t* sparsity,
                                          int num_local_block_rows,
                                          int num_remote_block_rows,
                                          int* block_sizes,
                                          ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);
  adj_graph_t* block_sparsity = adj_graph_new_with_block_sizes(sparsity, block_sizes);
  local_matrix_t* P = ilu_sparse_local_matrix_new(block_sparsity, ilu_params);
  adj_graph_free(block_sparsity);
  return var_cpr_pc_from_function_and_matrix(comm, F, NULL, context, dtor, 
                                             sparsity, num_local_block_rows, 
                                             num_remote_block_rows, block_sizes, P);
}
 
newton_pc_t* block_jacobi_cpr_pc_from_dae_function(MPI_Comm comm,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                   void* context,
                                                   void (*dtor)(void* context),
                                                   adj_graph_t* sparsity,
                                                   int num_local_block_rows,
                                                   int num_remote_block_rows,
                                                   int block_size)
{
  local_matrix_t* P = block_diagonal_matrix_new(num_local_block_rows, block_size);
  return cpr_pc_from_function_and_matrix(comm, NULL, F, context, dtor, 
                                         sparsity, num_local_block_rows, 
                                         num_remote_block_rows, block_size, P);
}

newton_pc_t* var_block_jacobi_cpr_pc_from_dae_function(MPI_Comm comm,
                                                       int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                       void* context,
                                                       void (*dtor)(void* context),
                                                       adj_graph_t* sparsity,
                                                       int num_local_block_rows,
                                                       int num_remote_block_rows,
                                                       int* block_sizes)
{
  local_matrix_t* P = var_block_diagonal_matrix_new(num_local_block_rows, block_sizes);
  return var_cpr_pc_from_function_and_matrix(comm, NULL, F, context, dtor, 
                                             sparsity, num_local_block_rows, 
                                             num_remote_block_rows, block_sizes, P);
}

newton_pc_t* lu_cpr_pc_from_dae_function(MPI_Comm comm,
                                         int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                         void* context,
                                         void (*dtor)(void* context),
                                         adj_graph_t* sparsity,
                                         int num_local_block_rows,
                                         int num_remote_block_rows,
                                         int block_size)
{
  adj_graph_t* block_sparsity = adj_graph_new_with_block_size(sparsity, block_size);
  local_matrix_t* P = sparse_local_matrix_new(block_sparsity);
  adj_graph_free(block_sparsity);
  return cpr_pc_from_function_and_matrix(comm, NULL, F, context, dtor, 
                                         sparsity, num_local_block_rows, 
                                         num_remote_block_rows, block_size, P);
}

newton_pc_t* var_lu_cpr_pc_from_dae_function(MPI_Comm comm,
                                             int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                             void* context,
                                             void (*dtor)(void* context),
                                             adj_graph_t* sparsity,
                                             int num_local_block_rows,
                                             int num_remote_block_rows,
                                             int* block_sizes)
{
  adj_graph_t* block_sparsity = adj_graph_new_with_block_sizes(sparsity, block_sizes);
  local_matrix_t* P = sparse_local_matrix_new(block_sparsity);
  adj_graph_free(block_sparsity);
  return var_cpr_pc_from_function_and_matrix(comm, NULL, F, context, dtor, 
                                             sparsity, num_local_block_rows, 
                                             num_remote_block_rows, block_sizes, P);
}

newton_pc_t* ilu_cpr_pc_from_dae_function(MPI_Comm comm,
                                          int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                          void* context,
                                          void (*dtor)(void* context),
                                          adj_graph_t* sparsity,
                                          int num_local_block_rows,
                                          int num_remote_block_rows,
                                          int block_size,
                                          ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);
  adj_graph_t* block_sparsity = adj_graph_new_with_block_size(sparsity, block_size);
  local_matrix_t* P = ilu_sparse_local_matrix_new(block_sparsity, ilu_params);
  adj_graph_free(block_sparsity);
  return cpr_pc_from_function_and_matrix(comm, NULL, F, context, dtor, 
                                         sparsity, num_local_block_rows, 
                                         num_remote_block_rows, block_size, P);
}

newton_pc_t* var_ilu_cpr_pc_from_dae_function(MPI_Comm comm,
                                              int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                              void* context,
                                              void (*dtor)(void* context),
                                              adj_graph_t* sparsity,
                                              int num_local_block_rows,
                                              int num_remote_block_rows,
                                              int* block_sizes,
                                              ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);
  adj_graph_t* block_sparsity = adj_graph_new_with_block_sizes(sparsity, block_sizes);
  local_matrix_t* P = ilu_sparse_local_matrix_new(block_sparsity, ilu_params);
  adj_graph_free(block_sparsity);
  return var_cpr_pc_from_function_and_matrix(comm, NULL, F, context, dtor, 
                                             sparsity, num_local_block_rows, 
                                             num_remote_block_rows, block_sizes, P);
}

