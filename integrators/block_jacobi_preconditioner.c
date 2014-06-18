// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/linear_algebra.h"
#include "core/sundials_helpers.h"
#include "integrators/block_jacobi_preconditioner.h"

// Block-diagonal matrix.
typedef struct
{
  int num_block_rows;
  int block_size;
  real_t* coeffs;
} bd_mat_t;

typedef struct 
{
  adj_graph_t* sparsity;
  adj_graph_coloring_t* coloring;
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*dae_F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* Fval);
  void* context;

  int num_block_rows;
  int block_size; 

  // Preconditioner matrix.
  bd_mat_t* P;

  // Work vectors.
  int num_work_vectors;
  real_t** work;

} block_jacobi_preconditioner_t;

static void bd_dtor(void* context)
{
  bd_mat_t* P = context;
  polymec_free(P->coeffs);
  polymec_free(P);
}

static bd_mat_t* block_jacobi_matrix(block_jacobi_preconditioner_t* precond)
{
  int bs = precond->block_size;
  bd_mat_t* P = polymec_malloc(sizeof(bd_mat_t));
  P->num_block_rows = precond->num_block_rows;
  P->block_size = bs;
  P->coeffs = polymec_malloc(sizeof(real_t) * P->num_block_rows * bs * bs);
  memset(P->coeffs, 0, sizeof(real_t) * P->num_block_rows * bs * bs);
  return P;
}

// Here's our finite difference implementation of the dF/dx matrix-vector 
// product. 
static void finite_diff_dFdx_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                               void* context, 
                               real_t t, 
                               real_t* x, 
                               real_t* xdot, 
                               int num_rows,
                               real_t* v, 
                               real_t** work, 
                               real_t* dFdx_v)
{
  real_t eps = sqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(t, x, xdot).
  // work[2] == x + eps*v
  // work[3] == F(t, x + eps*v)

  // x + eps*v -> work[2].
  for (int i = 0; i < num_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v, xdot) -> work[3].
  F(context, t, work[2], xdot, work[3]);

  // (F(t, x + eps*v, xdot) - F(t, x, xdot)) / eps -> (dF/dx) * v
  for (int i = 0; i < num_rows; ++i)
    dFdx_v[i] = (work[3][i] - work[1][i]) / eps;
}

// Here's the same finite difference calculation for dF/d(xdot).
static void finite_diff_dFdxdot_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                                  void* context, 
                                  real_t t, 
                                  real_t* x, 
                                  real_t* xdot, 
                                  int num_rows,
                                  real_t* v, 
                                  real_t** work, 
                                  real_t* dFdxdot_v)
{
  real_t eps = sqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(t, x, xdot).
  // work[2] == xdot + eps*v
  // work[3] == F(t, x, xdot + eps*v)

  // xdot + eps*v -> work[2].
  for (int i = 0; i < num_rows; ++i)
    work[2][i] = xdot[i] + eps*v[i];

  // F(t, x, xdot + eps*v) -> work[3].
  F(context, t, x, work[2], work[3]);

  // (F(t, x, xdot + eps*v) - F(t, x, xdot)) / eps -> (dF/dx) * v
  for (int i = 0; i < num_rows; ++i)
    dFdxdot_v[i] = (work[3][i] - work[1][i]) / eps;
}

static void add_Jv_into_bd_mat(adj_graph_t* graph,
                               adj_graph_coloring_t* coloring, 
                               int color, 
                               real_t factor,
                               real_t* Jv, 
                               bd_mat_t* J)
{
  int bs = J->block_size;
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    int block_col = i / bs;
    int c = i % bs;
    for (int j = block_col*bs; j < (block_col+1)*bs; ++j)
    {
      int r = j % bs;
      J->coeffs[block_col*bs*bs + c*bs + r] += factor * Jv[j];
    }
  }
}

// This function adapts non-DAE functions F(t, x) to DAE ones F(t, x, xdot).
static int F_adaptor(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval)
{
  ASSERT(xdot == NULL);

  // We are passed the actual preconditioner as our context pointer, so get the 
  // "real" one here.
  block_jacobi_preconditioner_t* pc = context;
  void* ctx = pc->context;
  return pc->F(ctx, t, x, Fval);
}

static void block_jacobi_preconditioner_setup(void* context, real_t alpha, real_t beta, real_t gamma,
                                              real_t t, real_t* x, real_t* xdot)
{
  block_jacobi_preconditioner_t* precond = context;
  adj_graph_t* graph = precond->sparsity;
  adj_graph_coloring_t* coloring = precond->coloring;
  real_t** work = precond->work;

  // Zero our preconditioner matrix and set the diagonal term.
  {
    bd_mat_t* P = precond->P;
    int n = P->num_block_rows;
    int bs = P->block_size;
    memset(P->coeffs, 0, sizeof(real_t) * n * bs * bs);
    for (int i = 0; i < n; ++i)
    {
      real_t* A = &P->coeffs[i*bs*bs];
      for (int j = 0; j < bs; ++j)
        A[j*bs+j] = alpha;
    }
  }

  int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* F_context;
  if (precond->dae_F != NULL)
  {
    ASSERT(xdot != NULL);
    F = precond->dae_F;
    F_context = precond->context;
  }
  else
  {
    ASSERT(xdot == NULL);
    F = F_adaptor;
    F_context = precond;
  }

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int num_rows = adj_graph_num_vertices(graph);
  ASSERT(num_rows == precond->block_size*precond->num_block_rows);
  real_t* Jv = polymec_malloc(sizeof(real_t) * num_rows);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * num_rows);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      work[0][i] = 1.0;

    // We evaluate F(t, x, xdot) and place it into work[1].
    F(F_context, t, x, xdot, work[1]);

    // Evaluate dF/dx and stash it in P.
    memset(Jv, 0, sizeof(real_t) * num_rows);
    finite_diff_dFdx_v(F, F_context, t, x, xdot, num_rows, work[0], work, Jv);
    add_Jv_into_bd_mat(graph, coloring, c, beta, Jv, precond->P);

    if (xdot != NULL)
    {
      // Now evaluate dF/d(xdot) and do the same.
      memset(Jv, 0, sizeof(real_t) * num_rows);
      finite_diff_dFdxdot_v(F, F_context, t, x, xdot, num_rows, work[0], work, Jv);
      add_Jv_into_bd_mat(graph, coloring, c, gamma, Jv, precond->P);
    }
  }
  polymec_free(Jv);
}


static bool block_jacobi_preconditioner_solve(void* context, real_t* B)
{
  block_jacobi_preconditioner_t* precond = context;
  int bs = precond->block_size;
  bd_mat_t* P = precond->P;

  bool success = false;
  for (int i = 0; i < precond->num_block_rows; ++i)
  {
    // Copy the block for this row into place.
    real_t Aij[bs*bs], bi[bs];
    memcpy(Aij, &P->coeffs[i*bs*bs], sizeof(real_t)*bs*bs);
    memcpy(bi, &B[i*bs], sizeof(real_t)*bs);

    // Replace each zero on the diagonal of Aij with a small number.
    static const real_t epsilon = 1e-25;
    for (int j = 0; j < bs; ++j)
    {
      if (Aij[bs*j+j] == 0.0)
        Aij[bs*j+j] = epsilon;
    }

    // Solve the linear system.
    int one = 1, ipiv[bs], info;
    dgesv(&bs, &one, Aij, &bs, ipiv, bi, &bs, &info);
    success = (info == 0);

    if (success)
    {
      // Copy the solution into place.
      memcpy(&B[i*bs], bi, sizeof(real_t)*bs);
    }
    else
    {
      ASSERT(info > 0);
      log_debug("block_jacobi_preconditioner_solve: call to dgesv failed for block row %d.", i);
      log_debug("(U is singular).", i);
      break;
    }
  }

  return success;
}

static void block_jacobi_preconditioner_fprintf(void* context, FILE* stream)
{
  block_jacobi_preconditioner_t* precond = context;
  bd_mat_t* P = precond->P;
  int n = P->num_block_rows;
  int bs = P->block_size;
  fprintf(stream, "Block Jacobian preconditioner matrix: (%d rows):", n * bs);
  for (int i = 0; i < n; ++i)
  {
    fprintf(stream, "\nRows %6d - %6d: ", i*bs, (i+1)*bs - 1);
    matrix_fprintf(&P->coeffs[i * bs * bs], bs, bs, stream);
  }
  fprintf(stream, "\n");
}

static void block_jacobi_preconditioner_dtor(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  for (int i = 0; i < precond->num_work_vectors; ++i)
    polymec_free(precond->work[i]);
  polymec_free(precond->work);
  adj_graph_coloring_free(precond->coloring);
  adj_graph_free(precond->sparsity);
  polymec_free(precond);
}

static preconditioner_t* general_block_jacobi_preconditioner_new(void* context,
                                                                 adj_graph_t* sparsity,
                                                                 int num_block_rows,
                                                                 int block_size)
{
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);

  block_jacobi_preconditioner_t* precond = polymec_malloc(sizeof(block_jacobi_preconditioner_t));

  // Do we have a block graph?
  int num_rows = adj_graph_num_vertices(sparsity);
  ASSERT((num_rows == num_block_rows) || (num_rows = block_size*num_block_rows));
  if (num_rows == num_block_rows)
    precond->sparsity = adj_graph_new_with_block_size(block_size, sparsity);
  else
    precond->sparsity = adj_graph_clone(sparsity);

  precond->coloring = adj_graph_coloring_new(precond->sparsity, SMALLEST_LAST);
  log_debug("Block Jacobi preconditioner: graph coloring produced %d colors.", 
            adj_graph_coloring_num_colors(precond->coloring));
  precond->context = context;
  precond->num_block_rows = num_block_rows;
  precond->block_size = block_size;
  precond->P = block_jacobi_matrix(precond);
  precond->F = NULL;
  precond->dae_F = NULL;

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = polymec_malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = polymec_malloc(sizeof(real_t) * precond->num_block_rows * precond->block_size);

  preconditioner_vtable vtable = {.setup = block_jacobi_preconditioner_setup,
                                  .solve = block_jacobi_preconditioner_solve,
                                  .fprintf = block_jacobi_preconditioner_fprintf,
                                  .dtor = block_jacobi_preconditioner_dtor};
  return preconditioner_new("Block Jacobi preconditioner", precond, vtable);
}

preconditioner_t* block_jacobi_preconditioner_new(void* context,
                                                  int (*F)(void* context, real_t t, real_t* x, real_t* func),
                                                  adj_graph_t* sparsity,
                                                  int num_block_rows,
                                                  int block_size)
{
  ASSERT(F != NULL);
  preconditioner_t* pc = general_block_jacobi_preconditioner_new(context, sparsity, num_block_rows, block_size);
  block_jacobi_preconditioner_t* bj = preconditioner_context(pc);
  bj->F = F;
  return pc;
}

preconditioner_t* block_jacobi_dae_preconditioner_new(void* context,
                                                      int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* func),
                                                      adj_graph_t* sparsity,
                                                      int num_block_rows,
                                                      int block_size)
{
  ASSERT(F != NULL);
  preconditioner_t* pc = general_block_jacobi_preconditioner_new(context, sparsity, num_block_rows, block_size);
  block_jacobi_preconditioner_t* bj = preconditioner_context(pc);
  bj->dae_F = F;
  return pc;
}

