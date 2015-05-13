// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h" // For UNIT_ROUNDOFF
#include "core/timer.h"
#include "integrators/cpr_differencer.h"

struct cpr_differencer_t
{
  MPI_Comm comm;

  // Function information.
  void* F_context;
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*F_dae)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void (*F_dtor)(void* context);

  // Matrix information.
  adj_graph_t* sparsity;
  int num_local_rows, num_remote_rows;
  adj_graph_coloring_t* coloring;
  int max_colors;

  // Work space.
  real_t* Jv;
  int num_work_vectors;
  real_t** work;
};

cpr_differencer_t* cpr_differencer_new(MPI_Comm comm,
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

  // Steal the graph.
  diff->sparsity = sparsity;

  // The number of local rows is known at this point.
  diff->num_local_rows = num_local_rows;

  // However, we can't know the actual number of remote block rows without 
  // knowing the specific communication pattern! However, it is safe to just 
  // multiply the given number of remote block rows by the maximum block size, 
  // since only the underlying RHS function will actually use the data.
  diff->num_remote_rows = num_remote_rows;

  // Assemble a graph coloring for any matrix we treat.
  diff->coloring = adj_graph_coloring_new(diff->sparsity, SMALLEST_LAST);

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

void cpr_differencer_free(cpr_differencer_t* diff)
{
  if ((diff->F_dtor != NULL) && (diff->F_context != NULL))
    diff->F_dtor(diff->F_context);
  adj_graph_coloring_free(diff->coloring);
  adj_graph_free(diff->sparsity);
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
  real_t eps = sqrt(UNIT_ROUNDOFF);
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

  // (F(t, x + eps*v, xdot) - F(t, x, xdot)) / eps -> (dF/dx) * v
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
  real_t eps = sqrt(UNIT_ROUNDOFF);
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

  // (F(t, x, xdot + eps*v) - F(t, x, xdot)) / eps -> (dF/dx) * v
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

void cpr_differencer_compute(cpr_differencer_t* diff, 
                             real_t alpha, real_t beta, real_t gamma,  
                             real_t t, real_t* x, real_t* xdot,
                             local_matrix_t* matrix)
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
    ASSERT(gamma == 0.0);
    ASSERT(xdot == NULL);
    F = F_adaptor;
    F_context = diff;
  }

  // First, zero the matrix.
  local_matrix_zero(matrix);

  // If all the coefficients are zero, we're finished!
  if ((alpha == 0.0) && (beta == 0.0) && (gamma == 0.0))
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Then set up an identity matrix.
  local_matrix_add_identity(matrix, alpha);

  // If beta and gamma are zero, we're finished!
  if ((beta == 0.0) && (gamma == 0.0))
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Keep track of the number of function evaluations.
  int num_F_evals = 0;

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  if ((alpha != 0.0) && (beta != 0.0) && (gamma != 0.0))
    log_debug("cpr_differencer: approximating J = %g * I + %g * dF/dx + %g * dF/d(xdot)...", alpha, beta, gamma);
  else if ((alpha == 0.0) && (beta != 0.0) && (gamma != 0.0))
    log_debug("cpr_differencer: approximating J = %g * dF/dx + %g * dF/d(xdot)...", beta, gamma);
  else if ((alpha == 0.0) && (beta == 0.0) && (gamma != 0.0))
    log_debug("cpr_differencer: approximating J = %g * dF/d(xdot)...", gamma);
  else if ((alpha != 0.0) && (beta != 0.0))
    log_debug("cpr_differencer: approximating J = %g * I + %g * dF/dx...", alpha, beta);
  else if ((alpha == 0.0) && (beta != 0.0))
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

    // Add the column vector J*v into our matrix.
    pos = 0;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      local_matrix_add_column_vector(matrix, beta, i, diff->Jv);

    if ((gamma != 0.0) && (xdot != NULL))
    {
      // Now evaluate dF/d(xdot) * d.
      memset(diff->Jv, 0, sizeof(real_t) * diff->num_local_rows);
      cpr_finite_diff_dFdxdot_v(F_context, F, t, x, xdot, diff->num_local_rows, 
                                diff->num_remote_rows, work[0], work, diff->Jv);
      ++num_F_evals;

      // Add in the column vector.
      pos = 0;
      while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
        local_matrix_add_column_vector(matrix, gamma, i, diff->Jv);
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
    if ((gamma != 0.0) && (xdot != NULL))
    {
      cpr_finite_diff_dFdxdot_v(F_context, F, t, x, xdot, diff->num_local_rows, 
                                diff->num_remote_rows, work[0], work, diff->Jv);
      ++num_F_evals;
    }
  }

  log_debug("cpr_differencer: Evaluated F %d times.", num_F_evals);
  STOP_FUNCTION_TIMER();
}

