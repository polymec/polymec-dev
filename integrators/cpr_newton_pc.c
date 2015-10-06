// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/block_diagonal_matrix.h"
#include "core/sparse_local_matrix.h"
#include "integrators/cpr_newton_pc.h"
#include "integrators/cpr_differencer.h"

// This preconditioner consists of a differencer and a given matrix.
typedef struct 
{
  cpr_differencer_t* diff;
  local_matrix_t* P;
} cpr_newton_pc_t;

static void cpr_newton_pc_compute_p(void* context, 
                                    real_t alpha, real_t beta, real_t gamma, 
                                    real_t t, real_t* x, real_t* xdot)
{
  cpr_newton_pc_t* pc = context;
  cpr_differencer_compute(pc->diff, alpha, beta, gamma, t, x, xdot, pc->P);
}

static bool cpr_newton_pc_solve(void* context, 
                                real_t t, real_t* x, real_t* xdot, real_t tolerance,
                                real_t* r, real_t* z, real_t* error_L2_norm)
{
  cpr_newton_pc_t* pc = context;
  bool solved = local_matrix_solve(pc->P, r, z);
  if (solved)
  {
    // Compute the L2 norm and measure against tolerance.
    int N = local_matrix_num_rows(pc->P);
    real_t Pz[N];
    local_matrix_matvec(pc->P, z, Pz);
    *error_L2_norm = 0.0;
    for (int i = 0; i < N; ++i)
      *error_L2_norm += (Pz[i]-r[i])*(Pz[i]-r[i]);
    *error_L2_norm = sqrt(*error_L2_norm);
    if (*error_L2_norm >= tolerance)
      solved = false;
  }
  return solved;
}

static void cpr_newton_pc_dtor(void* context)
{
  cpr_newton_pc_t* pc = context;
  cpr_differencer_free(pc->diff);
  local_matrix_free(pc->P);
  polymec_free(pc);
}

static newton_pc_t* cpr_newton_pc_from_function_and_matrix(const char* name, 
                                                           MPI_Comm comm,
                                                           void* context,
                                                           int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                           int (*dae_F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                           void (*dtor)(void* context),
                                                           newton_pc_side_t side,
                                                           adj_graph_t* sparsity,
                                                           int num_local_rows,
                                                           int num_remote_rows,
                                                           local_matrix_t* P)
{
  cpr_newton_pc_t* pc = polymec_malloc(sizeof(cpr_newton_pc_t));

  // DAE systems support only left preconditioning.
  ASSERT((side == NEWTON_PC_LEFT) || (dae_F == NULL));

  // Create a copy of the sparsity graph so the differencer can eat it.
  adj_graph_t* my_sparsity = adj_graph_clone(sparsity);

  pc->diff = cpr_differencer_new(comm, context, F, dae_F, dtor,
                                 my_sparsity, num_local_rows,
                                 num_remote_rows);
  pc->P = P;
  newton_pc_vtable vtable = {.compute_p = cpr_newton_pc_compute_p,
                             .solve = cpr_newton_pc_solve,
                             .dtor = cpr_newton_pc_dtor};
  return newton_pc_new(name, pc, vtable, side);
}

newton_pc_t* block_jacobi_cpr_newton_pc_from_function(MPI_Comm comm,
                                                      void* context,
                                                      int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                      void (*dtor)(void* context),
                                                      newton_pc_side_t side,
                                                      adj_graph_t* sparsity,
                                                      int num_local_block_rows,
                                                      int num_remote_block_rows,
                                                      int block_size)
{
  int num_local_rows = adj_graph_num_vertices(sparsity);
  ASSERT(num_local_rows == block_size * num_local_block_rows);
  int num_remote_rows = block_size * num_remote_block_rows;
  local_matrix_t* P = block_diagonal_matrix_new(num_local_block_rows, block_size);
  return cpr_newton_pc_from_function_and_matrix("Block-Jacobi CPR preconditioner", 
                                                comm, context, F, NULL, dtor, 
                                                side, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}
                                        
newton_pc_t* var_block_jacobi_cpr_newton_pc_from_function(MPI_Comm comm,
                                                          void* context,
                                                          int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                          void (*dtor)(void* context),
                                                          newton_pc_side_t side,
                                                          adj_graph_t* sparsity,
                                                          int num_local_block_rows,
                                                          int num_remote_block_rows,
                                                          int* block_sizes)
{
  int num_local_rows = adj_graph_num_vertices(sparsity);
  int alleged_num_local_rows = 0;
  int max_block_size = 1;
  for (int r = 0; r < num_local_block_rows; ++r)
  {
    ASSERT(block_sizes[r] >= 1);
    max_block_size = MAX(block_sizes[r], max_block_size);
    alleged_num_local_rows += block_sizes[r];
  }
  ASSERT(num_local_rows == alleged_num_local_rows);
  int num_remote_rows = max_block_size * num_remote_block_rows;
  local_matrix_t* P = var_block_diagonal_matrix_new(num_local_block_rows, block_sizes);
  return cpr_newton_pc_from_function_and_matrix("Block-Jacobi (var size) CPR preconditioner", 
                                                comm, context, F, NULL, dtor, 
                                                side, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}
 
newton_pc_t* lu_cpr_newton_pc_from_function(MPI_Comm comm,
                                            void* context,
                                            int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                            void (*dtor)(void* context),
                                            newton_pc_side_t side,
                                            adj_graph_t* sparsity,
                                            int num_local_rows,
                                            int num_remote_rows)
{
  ASSERT(num_local_rows == adj_graph_num_vertices(sparsity));
  local_matrix_t* P = sparse_local_matrix_new(sparsity);
  return cpr_newton_pc_from_function_and_matrix("LU CPR preconditioner", 
                                                comm, context, F, NULL, dtor, 
                                                side, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}
                                        
newton_pc_t* ilu_cpr_newton_pc_from_function(MPI_Comm comm,
                                             void* context,
                                             int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                             void (*dtor)(void* context),
                                             newton_pc_side_t side,
                                             adj_graph_t* sparsity,
                                             int num_local_rows,
                                             int num_remote_rows,
                                             ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);
  ASSERT(num_local_rows == adj_graph_num_vertices(sparsity));
  local_matrix_t* P = sparse_local_matrix_new(sparsity);
  sparse_local_matrix_use_ilu(P, ilu_params);
  return cpr_newton_pc_from_function_and_matrix("ILU CPR preconditioner", 
                                                comm, context, F, NULL, dtor, 
                                                side, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}
                                        
newton_pc_t* block_jacobi_cpr_newton_pc_from_dae_function(MPI_Comm comm,
                                                          void* context,
                                                          int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                          void (*dtor)(void* context),
                                                          adj_graph_t* sparsity,
                                                          int num_local_block_rows,
                                                          int num_remote_block_rows,
                                                          int block_size)
{
  int num_local_rows = adj_graph_num_vertices(sparsity);
  ASSERT(num_local_rows == block_size * num_local_block_rows);
  int num_remote_rows = block_size * num_remote_block_rows;
  local_matrix_t* P = block_diagonal_matrix_new(num_local_block_rows, block_size);
  return cpr_newton_pc_from_function_and_matrix("Block-Jacobi DAE CPR preconditioner", 
                                                comm, context, NULL, F, dtor, 
                                                NEWTON_PC_LEFT, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}

newton_pc_t* var_block_jacobi_cpr_newton_pc_from_dae_function(MPI_Comm comm,
                                                              void* context,
                                                              int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                              void (*dtor)(void* context),
                                                              adj_graph_t* sparsity,
                                                              int num_local_block_rows,
                                                              int num_remote_block_rows,
                                                              int* block_sizes)
{
  int num_local_rows = adj_graph_num_vertices(sparsity);
  int alleged_num_local_rows = 0;
  int max_block_size = 1;
  for (int r = 0; r < num_local_block_rows; ++r)
  {
    ASSERT(block_sizes[r] >= 1);
    max_block_size = MAX(block_sizes[r], max_block_size);
    alleged_num_local_rows += block_sizes[r];
  }
  ASSERT(num_local_rows == alleged_num_local_rows);
  int num_remote_rows = max_block_size * num_remote_block_rows;
  local_matrix_t* P = var_block_diagonal_matrix_new(num_local_block_rows, block_sizes);
  return cpr_newton_pc_from_function_and_matrix("Block-Jacobi (var size) DAE CPR preconditioner",
                                                comm, context, NULL, F, dtor, 
                                                NEWTON_PC_LEFT, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}

newton_pc_t* lu_cpr_newton_pc_from_dae_function(MPI_Comm comm,
                                                void* context,
                                                int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                void (*dtor)(void* context),
                                                adj_graph_t* sparsity,
                                                int num_local_rows,
                                                int num_remote_rows)
{
  ASSERT(num_local_rows == adj_graph_num_vertices(sparsity));
  local_matrix_t* P = sparse_local_matrix_new(sparsity);
  return cpr_newton_pc_from_function_and_matrix("LU DAE CPR preconditioner", 
                                                comm, context, NULL, F, dtor, 
                                                NEWTON_PC_LEFT, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}

newton_pc_t* ilu_cpr_newton_pc_from_dae_function(MPI_Comm comm,
                                                 void* context,
                                                 int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                 void (*dtor)(void* context),
                                                 adj_graph_t* sparsity,
                                                 int num_local_rows,
                                                 int num_remote_rows,
                                                 ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);
  ASSERT(num_local_rows == adj_graph_num_vertices(sparsity));
  local_matrix_t* P = sparse_local_matrix_new(sparsity);
  sparse_local_matrix_use_ilu(P, ilu_params);
  return cpr_newton_pc_from_function_and_matrix("ILU DAE CPR preconditioner", 
                                                comm, context, NULL, F, dtor, 
                                                NEWTON_PC_LEFT, sparsity, num_local_rows, 
                                                num_remote_rows, P);
}

bool newton_pc_is_cpr_newton_pc(newton_pc_t* cpr_newton_pc)
{
  // We just use the name.
  return (string_contains(newton_pc_name(cpr_newton_pc), "CPR preconditioner"));
}

local_matrix_t* cpr_newton_pc_matrix(newton_pc_t* cpr_newton_pc)
{
  // Must be a Curtis-Powell-Reed preconditioner!
  ASSERT(newton_pc_is_cpr_newton_pc(cpr_newton_pc));

  cpr_newton_pc_t* pc = newton_pc_context(cpr_newton_pc);
  return pc->P;
}
