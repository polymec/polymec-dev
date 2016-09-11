// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/bd_matrix.h"
#include "integrators/bj_newton_pc.h"

typedef struct 
{
  void* context;
  void (*compute_diag_block)(void* context, int i, 
                             real_t alpha, real_t beta, real_t gamma,
                             real_t t, real_t* x, real_t* x_dot, real_t* diag_block);
  void (*dtor)(void* context);
  int num_block_rows;
  bd_matrix_t* D;
} bj_pc_t;

static void bj_compute_p(void* context, 
                         real_t alpha, real_t beta, real_t gamma, 
                         real_t t, real_t* x, real_t* x_dot)
{
  bj_pc_t* pc = context;
  for (int i = 0; i < pc->num_block_rows; ++i)
  {
    real_t* J = bd_matrix_block(pc->D, i);
    pc->compute_diag_block(pc->context, i, alpha, beta, gamma, t, x, x_dot, J);
  }
}

static bool bj_solve(void* context, 
                     real_t t, real_t* x, real_t* xdot, real_t tolerance,
                     real_t* r, real_t* z, real_t* error_L2_norm)
{
  bj_pc_t* pc = context;
  bool solved = solve_bd_system(pc->D, r, z);
  if (solved)
  {
    // Compute the L2 norm and measure against tolerance.
    int N = bd_matrix_num_rows(pc->D);
    real_t Pz[N];
    bd_matrix_matvec(pc->D, z, Pz);
    *error_L2_norm = 0.0;
    for (int i = 0; i < N; ++i)
      *error_L2_norm += (Pz[i]-r[i])*(Pz[i]-r[i]);
    *error_L2_norm = sqrt(*error_L2_norm);
    if (*error_L2_norm >= tolerance)
      solved = false;
  }
  return solved;
}

static void bj_free(void* context)
{
  bj_pc_t* pc = context;
  bd_matrix_free(pc->D);
  polymec_free(pc);
}

newton_pc_t* bj_newton_pc_new(void* context,
                              void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                              void (*dtor)(void* context),
                              newton_pc_side_t side,
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
  pc->D = bd_matrix_new(num_block_rows, block_size);

  newton_pc_vtable vtable = {.compute_p = bj_compute_p,
                             .solve = bj_solve,
                             .dtor = bj_free};
  return newton_pc_new("Block Jacobi preconditioner", pc, vtable, side);

}
                                        
newton_pc_t* var_bj_newton_pc_new(void* context,
                                  void (*compute_diag_block)(void* context, int i, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* diag_block),
                                  void (*dtor)(void* context),
                                  newton_pc_side_t side,
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
  pc->D = var_bd_matrix_new(num_block_rows, block_sizes);

  newton_pc_vtable vtable = {.compute_p = bj_compute_p,
                             .solve = bj_solve,
                             .dtor = bj_free};
  return newton_pc_new("Variable Block Jacobi preconditioner", pc, vtable, side);

}
                                        
