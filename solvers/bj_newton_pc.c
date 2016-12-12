// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/bd_matrix.h"
#include "solvers/bj_newton_pc.h"

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

//------------------------------------------------------------------------
//            Curtis-Powell Reed block Jacobi preconditioners
//------------------------------------------------------------------------

// The cpr_differencer is an object that computes a Jacobian matrix for a 
// function F. The sparsity of the matrix is given by a graph, and the 
// differencer computes the matrix using the method of Curtis, Powell, and 
// Reed, in which the graph is colored so that the Jacobian can be computed 
// using finite differencing, with a minimum number of calls to F.
typedef struct 
{
  MPI_Comm comm;

  // Function information.
  void* F_context;
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*F_dae)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void (*F_dtor)(void* context);

  // Matrix information.
  int num_local_rows, num_remote_rows;
  adj_graph_coloring_t* coloring;
  int max_colors;

  // Work space.
  real_t* Jv;
  int num_work_vectors;
  real_t** work;
} cpr_differencer_t;

// Creates a new differencer.
static cpr_differencer_t* cpr_differencer_new(MPI_Comm comm,
                                              void* F_context,
                                              int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                              int (*F_dae)(void* context, real_t, real_t* x, real_t* xdot, real_t* Fval),
                                              void (*F_dtor)(void* context),
                                              adj_graph_t* sparsity,
                                              int num_local_rows,
                                              int num_remote_rows)
{
  ASSERT(num_local_rows > 0);
  ASSERT(num_remote_rows >= 0);

  // Exactly one of F and F_dae must be given.
  ASSERT((F != NULL) || (F_dae != NULL));
  ASSERT((F == NULL) || (F_dae == NULL));

  START_FUNCTION_TIMER();

  cpr_differencer_t* diff = polymec_malloc(sizeof(cpr_differencer_t));
  diff->comm = comm;
  diff->F_context = F_context;
  diff->F = F;
  diff->F_dae = F_dae;
  diff->F_dtor = F_dtor;

  // The number of local rows is known at this point.
  diff->num_local_rows = num_local_rows;

  // However, we can't know the actual number of remote block rows without 
  // knowing the specific communication pattern! However, it is safe to just 
  // multiply the given number of remote block rows by the maximum block size, 
  // since only the underlying RHS function will actually use the data.
  diff->num_remote_rows = num_remote_rows;

  // Assemble a graph coloring for any matrix we treat.
  diff->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);

  // Get the maximum number of colors on all MPI processes so that we can 
  // compute in lockstep.
  int num_colors = adj_graph_coloring_num_colors(diff->coloring);
  MPI_Allreduce(&num_colors, &diff->max_colors, 1, MPI_INT, MPI_MAX, diff->comm);
  log_debug("cpr_differencer: graph coloring produced %d colors.", diff->max_colors);

  // Make work vectors.
  diff->num_work_vectors = 4;
  diff->work = polymec_malloc(sizeof(real_t*) * diff->num_work_vectors);
  int N = diff->num_local_rows + diff->num_remote_rows;
  for (int i = 0; i < diff->num_work_vectors; ++i)
    diff->work[i] = polymec_malloc(sizeof(real_t) * N);
  diff->Jv = polymec_malloc(sizeof(real_t) * (diff->num_local_rows + diff->num_remote_rows));

  STOP_FUNCTION_TIMER();
  return diff;
}

static void cpr_differencer_free(cpr_differencer_t* diff)
{
  if ((diff->F_dtor != NULL) && (diff->F_context != NULL))
    diff->F_dtor(diff->F_context);
  for (int i = 0; i < diff->num_work_vectors; ++i)
    polymec_free(diff->work[i]);
  polymec_free(diff->work);
  adj_graph_coloring_free(diff->coloring);
  polymec_free(diff->Jv);
  polymec_free(diff);
}

// Here's our finite difference implementation of the dF/dx matrix-vector 
// product. 
static void cpr_finite_diff_dFdx_v(void* context, 
                                   int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                                   real_t t, 
                                   real_t* x, 
                                   real_t* xdot, 
                                   int num_local_rows,
                                   int num_remote_rows,
                                   real_t* v, 
                                   real_t** work, 
                                   real_t* dFdx_v)
{
  START_FUNCTION_TIMER();
  real_t eps = sqrt(REAL_EPSILON);
  real_t eps_inv = 1.0 / eps;

  // work[0] == v
  // work[1] contains F(t, x, xdot).
  // work[2] == x + eps*v
  // work[3] == F(t, x + eps*v)

  // x + eps*v -> work[2].
  for (int i = 0; i < num_local_rows + num_remote_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v, xdot) -> work[3].
  F(context, t, work[2], xdot, work[3]);

  // (F(t, x + eps*v, xdot) - F(t, x, xdot)) / eps -> dF/dx * v
  for (int i = 0; i < num_local_rows; ++i)
    dFdx_v[i] = (work[3][i] - work[1][i]) * eps_inv;
  STOP_FUNCTION_TIMER();
}

// Here's the same finite difference calculation for dF/d(xdot).
static void cpr_finite_diff_dFdxdot_v(void* context, 
                                      int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                                      real_t t, 
                                      real_t* x, 
                                      real_t* xdot, 
                                      int num_local_rows,
                                      int num_remote_rows,
                                      real_t* v, 
                                      real_t** work, 
                                      real_t* dFdxdot_v)
{
  START_FUNCTION_TIMER();
  real_t eps = sqrt(REAL_EPSILON);
  real_t eps_inv = 1.0 / eps;

  // work[0] == v
  // work[1] contains F(t, x, xdot).
  // work[2] == xdot + eps*v
  // work[3] == F(t, x, xdot + eps*v)

  // xdot + eps*v -> work[2].
  for (int i = 0; i < num_local_rows + num_remote_rows; ++i)
    work[2][i] = xdot[i] + eps*v[i];

  // F(t, x, xdot + eps*v) -> work[3].
  F(context, t, x, work[2], work[3]);

  // (F(t, x, xdot + eps*v) - F(t, x, xdot)) / eps -> dF/d(xdot) * v
  for (int i = 0; i < num_local_rows; ++i)
    dFdxdot_v[i] = (work[3][i] - work[1][i]) * eps_inv;
  STOP_FUNCTION_TIMER();
}

// This function adapts non-DAE functions F(t, x) to DAE ones F(t, x, xdot).
static int F_adaptor(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval)
{
  ASSERT(xdot == NULL);

  // We are passed the actual differencer as our context pointer, so get the 
  // "real" one here.
  cpr_differencer_t* diff = context;
  return diff->F(diff->F_context, t, x, Fval);
}

static void add_column_vector_to_matrix(void* matrix, real_t beta, int column, real_t* vec)
{
  bd_matrix_t* A = matrix;

  // We have to find the right block column.
  int block_col = 0, col_count = 0;
  while (true)
  {
    int block_size = bd_matrix_block_size(A, block_col);
    if ((col_count + block_size) > column) break;
    col_count += block_size;
    ++block_col;
  }

  int num_block_rows = bd_matrix_num_block_rows(A);
  if (block_col < num_block_rows)
  {
    int block_size = bd_matrix_block_size(A, block_col);
    int c = column % block_size;
    real_t* block = bd_matrix_block(A, block_col);
    for (int r = 0; r < block_size; ++r)
    {
      int j = col_count + r;
      block[c*block_size + r] += beta * vec[j];
    }
  }
}

static void cpr_differencer_compute(cpr_differencer_t* diff, 
                                    real_t alpha, real_t beta, real_t gamma,  
                                    real_t t, real_t* x, real_t* xdot,
                                    bd_matrix_t* matrix)
{
  START_FUNCTION_TIMER();
  adj_graph_coloring_t* coloring = diff->coloring;
  real_t** work = diff->work;

  // Normalize F.
  int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* F_context;
  if (diff->F_dae != NULL)
  {
    ASSERT(xdot != NULL);
    F = diff->F_dae;
    F_context = diff->F_context;
  }
  else
  {
    ASSERT(reals_equal(gamma, 0.0));
    ASSERT(xdot == NULL);
    F = F_adaptor;
    F_context = diff;
  }

  // First, zero the matrix.
  bd_matrix_zero(matrix);

  // If all the coefficients are zero, we're finished!
  if (reals_equal(alpha, 0.0) && reals_equal(beta, 0.0) && reals_equal(gamma, 0.0))
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Then set up an identity matrix.
  bd_matrix_add_identity(matrix, alpha);

  // If beta and gamma are zero, we're finished!
  if (reals_equal(beta, 0.0) && reals_equal(gamma, 0.0))
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Keep track of the number of function evaluations.
  int num_F_evals = 0;

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  if (!reals_equal(alpha, 0.0) && !reals_equal(beta, 0.0) && !reals_equal(gamma, 0.0))
    log_debug("cpr_differencer: approximating J = %g * I + %g * dF/dx + %g * dF/d(xdot)...", alpha, beta, gamma);
  else if (reals_equal(alpha, 0.0) && !reals_equal(beta, 0.0) && !reals_equal(gamma, 0.0))
    log_debug("cpr_differencer: approximating J = %g * dF/dx + %g * dF/d(xdot)...", beta, gamma);
  else if (reals_equal(alpha, 0.0) && reals_equal(beta, 0.0) && !reals_equal(gamma, 0.0))
    log_debug("cpr_differencer: approximating J = %g * dF/d(xdot)...", gamma);
  else if (!reals_equal(alpha, 0.0) && !reals_equal(beta, 0.0))
    log_debug("cpr_differencer: approximating J = %g * I + %g * dF/dx...", alpha, beta);
  else if (reals_equal(alpha, 0.0) && !reals_equal(beta, 0.0))
    log_debug("cpr_differencer: approximating J = %g * dF/dx...", beta);

  // Now iterate over all of the colors in our coloring. 
  int num_colors = adj_graph_coloring_num_colors(coloring);

  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * (diff->num_local_rows + diff->num_remote_rows));
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      work[0][i] = 1.0;

    // We evaluate F(t, x, xdot) and place it into work[1].
    F(F_context, t, x, xdot, work[1]);
    ++num_F_evals;

    // Evaluate dF/dx * d.
    memset(diff->Jv, 0, sizeof(real_t) * diff->num_local_rows);
    cpr_finite_diff_dFdx_v(F_context, F, t, x, xdot, diff->num_local_rows, 
                           diff->num_remote_rows, work[0], work, diff->Jv);
    ++num_F_evals;

    // Add the column vector dF/dx * d into our matrix.
    pos = 0;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      add_column_vector_to_matrix(matrix, beta, i, diff->Jv);

    if (!reals_equal(gamma, 0.0) && (xdot != NULL))
    {
      // Now evaluate dF/d(xdot) * d.
      memset(diff->Jv, 0, sizeof(real_t) * diff->num_local_rows);
      cpr_finite_diff_dFdxdot_v(F_context, F, t, x, xdot, diff->num_local_rows, 
                                diff->num_remote_rows, work[0], work, diff->Jv);
      ++num_F_evals;

      // Add in the column vector.
      pos = 0;
      while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
        add_column_vector_to_matrix(matrix, gamma, i, diff->Jv);
    }
  }

  // Now call the RHS functions in the same way as we would for all the colors
  // we don't have, up through the maximum number, so our neighboring 
  // processes can get exchanged data from us if they need it.
  int num_neighbor_colors = diff->max_colors - num_colors;
  for (int c = 0; c < num_neighbor_colors; ++c)
  {
    F(F_context, t, x, xdot, work[1]);
    ++num_F_evals;
    cpr_finite_diff_dFdx_v(F_context, F, t, x, xdot, diff->num_local_rows, 
                           diff->num_remote_rows, work[0], work, diff->Jv);
    ++num_F_evals;
    if (!reals_equal(gamma, 0.0) && (xdot != NULL))
    {
      cpr_finite_diff_dFdxdot_v(F_context, F, t, x, xdot, diff->num_local_rows, 
                                diff->num_remote_rows, work[0], work, diff->Jv);
      ++num_F_evals;
    }
  }

  log_debug("cpr_differencer: Evaluated F %d times.", num_F_evals);
  STOP_FUNCTION_TIMER();
}

// This preconditioner consists of our Curtis-Powell-Reed differencer and a given matrix.
typedef struct 
{
  cpr_differencer_t* diff;
  bd_matrix_t* P;
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
  bool solved = solve_bd_system(pc->P, r, z);
  if (solved)
  {
    // Compute the L2 norm and measure against tolerance.
    int N = bd_matrix_num_rows(pc->P);
    real_t Pz[N];
    bd_matrix_matvec(pc->P, z, Pz);
    *error_L2_norm = 0.0;
    polymec_suspend_fpe();
    for (int i = 0; i < N; ++i)
      *error_L2_norm += (Pz[i]-r[i])*(Pz[i]-r[i]);
    *error_L2_norm = sqrt(*error_L2_norm);
    polymec_restore_fpe();
    if (*error_L2_norm >= tolerance)
      solved = false;
  }
  return solved;
}

static void cpr_newton_pc_dtor(void* context)
{
  cpr_newton_pc_t* pc = context;
  cpr_differencer_free(pc->diff);
  bd_matrix_free(pc->P);
  polymec_free(pc);
}

static newton_pc_t* cpr_bj_newton_pc_from_function(MPI_Comm comm,
                                                   void* context,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                   int (*dae_F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                   void (*dtor)(void* context),
                                                   newton_pc_side_t side,
                                                   adj_graph_t* sparsity,
                                                   int num_local_rows,
                                                   int num_remote_rows,
                                                   bd_matrix_t* P)
{
  cpr_newton_pc_t* pc = polymec_malloc(sizeof(cpr_newton_pc_t));

  // DAE systems support only left preconditioning.
  ASSERT((side == NEWTON_PC_LEFT) || (dae_F == NULL));

  // Create a copy of the sparsity graph so the differencer can eat it.
  adj_graph_t* my_sparsity = adj_graph_clone(sparsity);
  pc->diff = cpr_differencer_new(comm, context, F, dae_F, dtor,
                                 my_sparsity, num_local_rows,
                                 num_remote_rows);
  adj_graph_free(my_sparsity);
  pc->P = P;
  newton_pc_vtable vtable = {.compute_p = cpr_newton_pc_compute_p,
                             .solve = cpr_newton_pc_solve,
                             .dtor = cpr_newton_pc_dtor};
  return newton_pc_new("Curtis-Powell-Reed block-Jacobi preconditioner", pc, vtable, side);
}

newton_pc_t* cpr_bj_newton_pc_new(MPI_Comm comm,
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
  bd_matrix_t* P = bd_matrix_new(num_local_block_rows, block_size);
  return cpr_bj_newton_pc_from_function(comm, context, F, NULL, dtor, 
                                        side, sparsity, num_local_rows, 
                                        num_remote_rows, P);
}
                                        
newton_pc_t* var_cpr_bj_newton_pc_new(MPI_Comm comm,
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
  bd_matrix_t* P = var_bd_matrix_new(num_local_block_rows, block_sizes);
  return cpr_bj_newton_pc_from_function(comm, context, F, NULL, dtor, 
                                        side, sparsity, num_local_rows, 
                                        num_remote_rows, P);
}

newton_pc_t* dae_cpr_bj_newton_pc_new(MPI_Comm comm,
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
  bd_matrix_t* P = bd_matrix_new(num_local_block_rows, block_size);
  return cpr_bj_newton_pc_from_function(comm, context, NULL, F, dtor, 
                                        NEWTON_PC_LEFT, sparsity, num_local_rows, 
                                        num_remote_rows, P);
}

newton_pc_t* var_dae_cpr_bj_newton_pc_new(MPI_Comm comm,
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
  bd_matrix_t* P = var_bd_matrix_new(num_local_block_rows, block_sizes);
  return cpr_bj_newton_pc_from_function(comm, context, NULL, F, dtor, 
                                        NEWTON_PC_LEFT, sparsity, num_local_rows, 
                                        num_remote_rows, P);
}
                                        
