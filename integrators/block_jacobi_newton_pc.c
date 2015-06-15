// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/block_diagonal_matrix.h"
#include "integrators/block_jacobi_newton_pc.h"

typedef struct 
{
  void* context;
  void (*compute_diag_block)(void* context, int i, 
                             real_t alpha, real_t beta, real_t gamma,
                             real_t t, real_t* x, real_t* x_dot, real_t* diag_block);
  void (*dtor)(void* context);
  int num_block_rows;
  local_matrix_t* D;
} bj_pc_t;

static void bj_compute_p(void* context, 
                         real_t alpha, real_t beta, real_t gamma, 
                         real_t t, real_t* x, real_t* x_dot)
{
  bj_pc_t* pc = context;
  for (int i = 0; i < pc->num_block_rows; ++i)
  {
    real_t* J = block_diagonal_matrix_row(pc->D, i);
    pc->compute_diag_block(pc->context, i, alpha, beta, gamma, t, x, x_dot, J);
  }
}

static bool bj_solve(void* context, 
                     real_t t, real_t* x, real_t* xdot,
                     real_t* r, real_t* z)
{
  bj_pc_t* pc = context;
  return local_matrix_solve(pc->D, r, z);
}

static void bj_free(void* context)
{
  bj_pc_t* pc = context;
  local_matrix_free(pc->D);
  polymec_free(pc);
}

newton_pc_t* block_jacobi_newton_pc_new(void* context,
                                        void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                                        void (*dtor)(void* context),
                                        int num_block_rows,
                                        int block_size)
{
  ASSERT(compute_diag_block != NULL);
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);

  bj_pc_t* pc = polymec_malloc(sizeof(bj_pc_t));
  pc->context = context;
  pc->compute_diag_block = compute_diag_block;
  pc->dtor = dtor;
  pc->num_block_rows = num_block_rows;
  pc->D = block_diagonal_matrix_new(num_block_rows, block_size);

  newton_pc_vtable vtable = {.compute_p = bj_compute_p,
                             .solve = bj_solve,
                             .dtor = bj_free};
  return newton_pc_new("Block Jacobi preconditioner", pc, vtable);

}
                                        
newton_pc_t* var_block_jacobi_newton_pc_new(void* context,
                                            void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                                            void (*dtor)(void* context),
                                            int num_block_rows,
                                            int* block_sizes)
{
  ASSERT(compute_diag_block != NULL);
  ASSERT(num_block_rows > 0);
  ASSERT(block_sizes != NULL);

  bj_pc_t* pc = polymec_malloc(sizeof(bj_pc_t));
  pc->context = context;
  pc->compute_diag_block = compute_diag_block;
  pc->dtor = dtor;
  pc->num_block_rows = num_block_rows;
  pc->D = var_block_diagonal_matrix_new(num_block_rows, block_sizes);

  newton_pc_vtable vtable = {.compute_p = bj_compute_p,
                             .solve = bj_solve,
                             .dtor = bj_free};
  return newton_pc_new("Variable Block Jacobi preconditioner", pc, vtable);

}
                                        
