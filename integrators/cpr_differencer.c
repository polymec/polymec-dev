// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "slu_ddefs.h"
#include "slu_util.h"
#include "core/sundials_helpers.h" // For UNIT_ROUNDOFF
#include "core/linear_algebra.h"
#include "core/array_utils.h"
#include "integrators/cpr_differencer.h"

struct cpr_differencer_t
{
  MPI_Comm comm;

  // Function information.
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*F_dae)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* F_context;

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
                                       int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                       int (*F_dae)(void* context, real_t, real_t* x, real_t* xdot, real_t* Fval),
                                       void* F_context,
                                       adj_graph_t* sparsity,
                                       int num_local_block_rows,
                                       int num_remote_block_rows,
                                       int block_size)
{
  ASSERT(num_local_block_rows > 0);
  ASSERT(block_size >= 1);

  int block_sizes[num_local_block_rows];
  for (int i = 0; i < num_local_block_rows; ++i)
    block_sizes[i] = block_size;

  return var_cpr_differencer_new(comm, F, F_dae, F_context, sparsity, num_local_block_rows,
                                 num_remote_block_rows, block_sizes);
}

cpr_differencer_t* var_cpr_differencer_new(MPI_Comm comm,
                                           int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                           int (*F_dae)(void* context, real_t, real_t* x, real_t* xdot, real_t* Fval),
                                           void* F_context,
                                           adj_graph_t* sparsity,
                                           int num_local_block_rows,
                                           int num_remote_block_rows,
                                           int* block_sizes)
{
  ASSERT(num_local_block_rows > 0);
  ASSERT(num_remote_block_rows >= 0);
  ASSERT(block_sizes != NULL);
#ifndef NDEBUG
  for (int i = 0; i < num_local_block_rows; ++i)
  {
    ASSERT(block_sizes[i] > 0);
  }
#endif

  // Exactly one of F and F_dae must be given.
  ASSERT((F != NULL) || (F_dae != NULL));
  ASSERT((F == NULL) || (F_dae == NULL));

  cpr_differencer_t* diff = polymec_malloc(sizeof(cpr_differencer_t));
  diff->comm = comm;
  diff->F_context = F_context;
  diff->F = F;
  diff->F_dae = F_dae;
  diff->Jv = polymec_malloc(sizeof(real_t) * (diff->num_local_rows + diff->num_remote_rows));

  // Do we have a block graph?
  int num_local_rows = adj_graph_num_vertices(sparsity);
  int alleged_num_local_rows = 0;
  int max_block_size = 1;
  for (int r = 0; r < num_local_block_rows; ++r)
  {
    ASSERT(block_sizes[r] >= 1);
    max_block_size = MAX(block_sizes[r], max_block_size);
    alleged_num_local_rows += block_sizes[r];
  }
  ASSERT((num_local_rows == num_local_block_rows) || (num_local_rows == alleged_num_local_rows)); 
  if (num_local_rows == num_local_block_rows)
  {
    // We were given the number of vertices in the graph as the number of 
    // block rows, so we create a graph with a block size of 1.
    diff->sparsity = adj_graph_new_with_block_sizes(block_sizes, sparsity);
    ASSERT(adj_graph_num_vertices(diff->sparsity) == alleged_num_local_rows);
  }
  else
  {
    // The number of vertices in the graph is the number of degrees of freedom
    // in the solution, so we don't need to create
    diff->sparsity = adj_graph_clone(sparsity);
  }

  // The number of local rows is known at this point.
  diff->num_local_rows = alleged_num_local_rows;

  // However, we can't know the actual number of remote block rows without 
  // knowing the specific communication pattern! However, it is safe to just 
  // multiply the given number of remote block rows by the maximum block size, 
  // since only the underlying RHS function will actually use the data.
  diff->num_remote_rows = num_remote_block_rows * max_block_size;

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

  return diff;
}

void cpr_differencer_free(cpr_differencer_t* diff)
{
  adj_graph_coloring_free(diff->coloring);
  adj_graph_free(diff->sparsity);
  polymec_free(diff->Jv);
  polymec_free(diff);
}

// Here's our finite difference implementation of the dF/dx matrix-vector 
// product. 
static void cpr_finite_diff_dFdx_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                                   void* context, 
                                   real_t t, 
                                   real_t* x, 
                                   real_t* xdot, 
                                   int num_local_rows,
                                   int num_remote_rows,
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
  for (int i = 0; i < num_local_rows + num_remote_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v, xdot) -> work[3].
  F(context, t, work[2], xdot, work[3]);

  // (F(t, x + eps*v, xdot) - F(t, x, xdot)) / eps -> (dF/dx) * v
  for (int i = 0; i < num_local_rows; ++i)
    dFdx_v[i] = (work[3][i] - work[1][i]) / eps;
}

// Here's the same finite difference calculation for dF/d(xdot).
static void cpr_finite_diff_dFdxdot_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                                      void* context, 
                                      real_t t, 
                                      real_t* x, 
                                      real_t* xdot, 
                                      int num_local_rows,
                                      int num_remote_rows,
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
  for (int i = 0; i < num_local_rows + num_remote_rows; ++i)
    work[2][i] = xdot[i] + eps*v[i];

  // F(t, x, xdot + eps*v) -> work[3].
  F(context, t, x, work[2], work[3]);

  // (F(t, x, xdot + eps*v) - F(t, x, xdot)) / eps -> (dF/dx) * v
  for (int i = 0; i < num_local_rows; ++i)
    dFdxdot_v[i] = (work[3][i] - work[1][i]) / eps;
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
    return;

  // Then set up an identity matrix.
  local_matrix_add_identity(matrix, alpha);

  // If beta and gamma are zero, we're finished!
  if ((beta == 0.0) && (gamma == 0.0))
    return;

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
    cpr_finite_diff_dFdx_v(F, F_context, t, x, xdot, diff->num_local_rows, 
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
      cpr_finite_diff_dFdxdot_v(F, F_context, t, x, xdot, diff->num_local_rows, 
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
    cpr_finite_diff_dFdx_v(F, F_context, t, x, xdot, diff->num_local_rows, 
                           diff->num_remote_rows, work[0], work, diff->Jv);
    ++num_F_evals;
    if ((gamma != 0.0) && (xdot != NULL))
    {
      cpr_finite_diff_dFdxdot_v(F, F_context, t, x, xdot, diff->num_local_rows, 
                                diff->num_remote_rows, work[0], work, diff->Jv);
      ++num_F_evals;
    }
  }

  log_debug("cpr_differencer: Evaluated F %d times.", num_F_evals);
}

struct local_matrix_t
{
  char* name;
  void* context;
  local_matrix_vtable vtable;
};

local_matrix_t* local_matrix_new(const char* name,
                                 void* context,
                                 local_matrix_vtable vtable)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_column_vector != NULL);
  ASSERT(vtable.add_identity != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(vtable.fprintf != NULL);
  local_matrix_t* matrix = polymec_malloc(sizeof(local_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  return matrix;
}
 
void local_matrix_free(local_matrix_t* matrix)
{
  if ((matrix->vtable.dtor != NULL) && (matrix->context != NULL))
    matrix->vtable.dtor(matrix->context);
  string_free(matrix->name);
  polymec_free(matrix);
}

char* local_matrix_name(local_matrix_t* matrix)
{
  return matrix->name;
}

void local_matrix_zero(local_matrix_t* matrix)
{
  matrix->vtable.zero(matrix->context);
}

void local_matrix_add_column_vector(local_matrix_t* matrix, 
                                    real_t scale_factor,
                                    int column,
                                    real_t* column_vector)
{
  matrix->vtable.add_column_vector(matrix->context,
                                   scale_factor,
                                   column,
                                   column_vector);
}

void local_matrix_add_identity(local_matrix_t* matrix, 
                               real_t scale_factor)
{
  matrix->vtable.add_identity(matrix->context, scale_factor);
}

bool local_matrix_solve(local_matrix_t* matrix, 
                       real_t* B,
                       real_t* x)
{
  return matrix->vtable.solve(matrix->context, B, x);
}

void local_matrix_fprintf(local_matrix_t* matrix, FILE* stream)
{
  matrix->vtable.fprintf(matrix->context, stream);
}

typedef struct 
{
  int num_block_rows;
  int *D_offsets, *B_offsets; // For variable block sizes.
  real_t* D;
  bool constant_block_size;
} bdm_t;

static void bdm_zero(void* context)
{
  bdm_t* A = context;
  memset(A->D, 0, sizeof(real_t) * A->D_offsets[A->num_block_rows]);
}

static void bdm_add_identity(void* context, real_t scale_factor)
{
  bdm_t* A = context;
  for (int i = 0; i < A->num_block_rows; ++i)
  {
    int bs = A->B_offsets[i+1] - A->B_offsets[i];
    int offset = A->D_offsets[i];
    for (int j = 0; j < bs; ++j)
      A->D[offset+bs*j+j] += scale_factor;
  }
}

static void bdm_add_column_vector(void* context,
                                  real_t scale_factor,
                                  int column,
                                  real_t* column_vector)
{
  bdm_t* A = context;
  real_t* D = A->D;
  int i = column;
  if (i < A->B_offsets[A->num_block_rows])
  {
    int bs = A->B_offsets[i+1] - A->B_offsets[i];
    int block_col = i / bs;
    int c = i % bs;
    for (int j = block_col*bs; j < (block_col+1)*bs; ++j)
    {
      int r = j % bs;
      D[A->D_offsets[block_col] + c*bs + r] += scale_factor * column_vector[j];
    }
  }
}

static bool bdm_solve(void* context, real_t* B, real_t* x)
{
  bdm_t* A = context;
  real_t* D = A->D;

  bool success = false;
  for (int i = 0; i < A->num_block_rows; ++i)
  {
    // Copy the block for this row into place.
    int bs = A->B_offsets[i+1] - A->B_offsets[i];
    int D_offset = A->D_offsets[i], B_offset = A->B_offsets[i];
    real_t Aij[bs*bs], bi[bs];
    memcpy(Aij, &D[D_offset], sizeof(real_t)*bs*bs);
    memcpy(bi, &B[B_offset], sizeof(real_t)*bs);

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
      memcpy(&x[B_offset], bi, sizeof(real_t)*bs);
    }
    else
    {
      ASSERT(info > 0);
      log_debug("bdm_solve: call to dgesv failed for block row %d.", i);
      log_debug("bdm_solve: (U is singular).", i);
      break;
    }

    D_offset += bs*bs;
    B_offset += bs;
  }

  return success;
}

static void bdm_fprintf(void* context, FILE* stream)
{
  bdm_t* A = context;
  int N = A->num_block_rows;
  real_t* D = A->D;
  fprintf(stream, "\nBlock diagonal matrix P:\n");
  for (int i = 0; i < N; ++i)
  {
    fprintf(stream, "%d: [", i);
    int bs = A->B_offsets[i+1] - A->B_offsets[i];
    int offset = A->D_offsets[i];
    for (int ii = 0; ii < bs; ++ii)
    {
      for (int jj = 0; jj < bs; ++jj)
        fprintf(stream, "%g ", D[offset+bs*ii+jj]);
      if (ii < (bs - 1))
        fprintf(stream, "; ");
    }
    fprintf(stream, "]\n");
  }
}

static void bdm_dtor(void* context)
{
  bdm_t* A = context;
  polymec_free(A->D);
  polymec_free(A->D_offsets);
  polymec_free(A->B_offsets);
  polymec_free(A);
}

local_matrix_t* block_diagonal_matrix_new(int num_block_rows,
                                          int block_size)
{
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);
  int block_sizes[num_block_rows];
  for (int i = 0; i < num_block_rows; ++i)
    block_sizes[i] = block_size;
  return var_block_diagonal_matrix_new(num_block_rows, block_sizes);
}

local_matrix_t* var_block_diagonal_matrix_new(int num_block_rows,
                                              int* block_sizes)
{
  bdm_t* A = polymec_malloc(sizeof(bdm_t));
  A->num_block_rows = num_block_rows;
  A->D_offsets = polymec_malloc(sizeof(int) * (num_block_rows+1));
  A->B_offsets = polymec_malloc(sizeof(int) * (num_block_rows+1));
  A->D_offsets[0] = A->B_offsets[0] = 0;
  A->constant_block_size = true;
  int bs0 = -1;
  for (int i = 0; i < num_block_rows; ++i)
  {
    int bs = block_sizes[i];
    if (bs0 == -1)
      bs0 = bs;
    else if (bs != bs0)
      A->constant_block_size = false;
    ASSERT(bs >= 1);
    A->D_offsets[i+1] = A->D_offsets[i] + bs*bs;
    A->B_offsets[i+1] = A->B_offsets[i] + bs;
  }
  int N = A->D_offsets[A->num_block_rows];
  A->D = polymec_malloc(sizeof(real_t) * N);

  char name[1024];
  if (A->constant_block_size)
    snprintf(name, 1024, "Block diagonal matrix (bs = %d)", bs0);
  else
    snprintf(name, 1024, "Variable block diagonal matrix");
  local_matrix_vtable vtable = {.dtor = bdm_dtor,
                                .zero = bdm_zero,
                                .add_identity = bdm_add_identity,
                                .add_column_vector = bdm_add_column_vector,
                                .solve = bdm_solve,
                                .fprintf = bdm_fprintf};
  return local_matrix_new(name, A, vtable);
}

typedef struct 
{
  int N; // Number of rows in matrix.

  // Personal copy of the sparsity graph.
  adj_graph_t* sparsity;

  // Matrix data.
  int* etree;
  SuperMatrix rhs, X, L, U;
  real_t *rhs_data, *X_data;
  int *rperm, *cperm;
  superlu_options_t options;
  SuperLUStat_t stat;
  SuperMatrix* A; // <-- matrix itself!.

  // ILU parameters (if any).
  ilu_params_t* ilu_params;
} slm_t;

static void supermatrix_free(SuperMatrix* matrix)
{
  switch(matrix->Stype) {
  case SLU_DN:
    Destroy_Dense_Matrix(matrix);
    break;
  case SLU_NR:
    Destroy_CompRow_Matrix(matrix);
    break;
  case SLU_NC:
    Destroy_CompCol_Matrix(matrix);
    break;
  default:
    printf("ERROR: unknown matrix type passed to supermatrix_free!");
    break;
  }
  polymec_free(matrix);
}

static SuperMatrix* supermatrix_new(adj_graph_t* graph)
{
  SuperMatrix* A = polymec_malloc(sizeof(SuperMatrix));
  
  // Fetch sparsity information from the graph.
  int* edge_offsets = adj_graph_edge_offsets(graph);
  int num_rows = adj_graph_num_vertices(graph);
  int num_edges = edge_offsets[num_rows];
  int num_nz = num_rows + num_edges;

  // Translate the sparsity information in the graph to column indices and 
  // row pointers. Since the graph is effectively stored in compressed 
  // row format, we will use SuperLU's compressed row matrix format.
  int* row_ptrs = intMalloc(num_rows + 1);
  int* col_indices = intMalloc(num_nz);
  int* edges = adj_graph_adjacency(graph);
  int offset = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    row_ptrs[i] = offset;
    col_indices[offset++] = i; // diagonal column

    // Off-diagonal columns.
    for (int j = edge_offsets[i]; j < edge_offsets[i+1]; ++j)
      col_indices[offset++] = edges[j];
    int_qsort(&col_indices[row_ptrs[i]+1], edge_offsets[i+1]-edge_offsets[i]);
  }
  row_ptrs[num_rows] = offset;
  ASSERT(offset == num_nz);

  // Create zeros for the matrix.
  real_t* mat_zeros = SUPERLU_MALLOC(sizeof(real_t) * num_nz);
  memset(mat_zeros, 0, sizeof(real_t) * num_nz);

  // Hand over these resources to create the Supermatrix.
  dCreate_CompCol_Matrix(A, num_rows, num_rows, num_nz, 
                         mat_zeros, col_indices, row_ptrs, 
                         SLU_NC, SLU_D, SLU_GE);
  return A;
}

static void slm_zero(void* context)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;
  int num_cols = A->ncol;
  NCformat* data = A->Store;
  real_t* Aij = data->nzval;
  memset(Aij, 0, data->colptr[num_cols]);
}

static void slm_add_identity(void* context, real_t scale_factor)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;
  int num_cols = A->ncol;
  NCformat* data = A->Store;
  real_t* Aij = data->nzval;

  // Add in the diagonal values.
  for (int j = 0; j < num_cols; ++j)
  {
    int col_index = data->colptr[j];
    Aij[col_index] += scale_factor;
  }
}

static void slm_add_column_vector(void* context, 
                                  real_t scale_factor,
                                  int column,
                                  real_t* column_vector)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;

  NCformat* data = A->Store;
  real_t* Jij = data->nzval;

  int i = column;

  if (i >= mat->N) return;

  // Add in the diagonal element.
  Jij[data->colptr[i]] += scale_factor * column_vector[i];
      
  // Add in off-diagonal column values.
  int pos = 0, j;
  while (adj_graph_next_edge(mat->sparsity, i, &pos, &j))
  {
    int col_index = data->colptr[i];
    size_t num_rows = data->colptr[i+1] - col_index;
    int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, j);
    ASSERT(entry != NULL);
    size_t offset = entry - &data->rowind[col_index];
    Jij[data->colptr[i] + offset] += scale_factor * column_vector[j];
  }
}

static bool slm_solve(void* context, real_t* B, real_t* x)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;

  // Copy B to the rhs vector.
  DNformat* rhs = mat->rhs.Store;
  memcpy(rhs->nzval, B, sizeof(real_t) * mat->N);

  if (mat->cperm == NULL)
  {
    mat->cperm = intMalloc(mat->N);
    mat->rperm = intMalloc(mat->N);
  }
  else if (mat->options.Fact != SamePattern)
  {
    // Jettison the existing factorization.
    Destroy_SuperNode_Matrix(&mat->L);
    Destroy_CompCol_Matrix(&mat->U);
  }

  // Do the solve.
  int info;
  dgssv(&mat->options, A, mat->cperm, mat->rperm, 
        &mat->L, &mat->U, &mat->rhs, &mat->stat, &info);
  // SuperLU claims that it can re-use factorization info with the same 
  // sparsity pattern, but this appears not to be the case at the moment.
  // mat->options.Fact = SamePattern;

  bool success = (info == 0);

  if (success)
  {
    // Copy the rhs vector to B.
    memcpy(B, rhs->nzval, sizeof(real_t) * mat->N);
  }

  return success;
}

static void slm_fprintf(void* context, FILE* stream)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;
  fprintf(stream, "\nCompCol matrix:\n");
  fprintf(stream, "Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
  int n = A->ncol;
  NCformat* Astore = (NCformat *) A->Store;
  real_t* dp = Astore->nzval;
  fprintf(stream, "nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
  fprintf(stream, "nzval: ");
  for (int i = 0; i < Astore->colptr[n]; ++i) 
    fprintf(stream, "%f  ", (double)dp[i]);
  fprintf(stream, "\nrowind: ");
  for (int i = 0; i < Astore->colptr[n]; ++i) 
    fprintf(stream, "%d  ", Astore->rowind[i]);
  fprintf(stream, "\ncolptr: ");
  for (int i = 0; i <= n; ++i) 
    fprintf(stream, "%d  ", Astore->colptr[i]);
  fprintf(stream, "\n");
}

static void slm_dtor(void* context)
{
  slm_t* mat = context;
  if (mat->cperm != NULL)
  {
    SUPERLU_FREE(mat->cperm);
    SUPERLU_FREE(mat->rperm);
    Destroy_SuperNode_Matrix(&mat->L);
    Destroy_CompCol_Matrix(&mat->U);
  }
  Destroy_SuperMatrix_Store(&mat->rhs);
  polymec_free(mat->rhs_data);
  Destroy_SuperMatrix_Store(&mat->X);
  supermatrix_free(mat->A);
  polymec_free(mat->X_data);
  StatFree(&mat->stat);
  if (mat->etree != NULL)
    polymec_free(mat->etree);
  mat->ilu_params = NULL;
  polymec_free(mat);
}

local_matrix_t* sparse_local_matrix_new(adj_graph_t* sparsity,
                                        int num_block_rows,
                                        int block_size)
{
  ASSERT(block_size > 0);
  int block_sizes[num_block_rows];
  for (int i = 0; i < num_block_rows; ++i)
    block_sizes[i] = block_size;
  return var_sparse_local_matrix_new(sparsity, num_block_rows, block_sizes);
}

local_matrix_t* var_sparse_local_matrix_new(adj_graph_t* sparsity,
                                            int num_block_rows,
                                            int* block_sizes)
{
  return var_ilu_sparse_local_matrix_new(sparsity, num_block_rows, block_sizes, NULL);
}

local_matrix_t* ilu_sparse_local_matrix_new(adj_graph_t* sparsity,
                                            int num_block_rows,
                                            int block_size,
                                            ilu_params_t* ilu_params)
{
  ASSERT(block_size > 0);
  int block_sizes[num_block_rows];
  for (int i = 0; i < num_block_rows; ++i)
    block_sizes[i] = block_size;
  return var_ilu_sparse_local_matrix_new(sparsity, num_block_rows, 
                                         block_sizes, ilu_params);
}

local_matrix_t* ilu_var_sparse_local_matrix_new(adj_graph_t* sparsity,
                                                int num_block_rows,
                                                int* block_sizes, 
                                                ilu_params_t* ilu_params)
{
  slm_t* mat = polymec_malloc(sizeof(slm_t));
  mat->sparsity = adj_graph_clone(sparsity); // MINE!
  mat->ilu_params = NULL;
  mat->A = supermatrix_new(sparsity);

  // Solver data.
  mat->N = adj_graph_num_vertices(sparsity);
  mat->rhs_data = polymec_malloc(sizeof(real_t) * mat->N);
  dCreate_Dense_Matrix(&mat->rhs, mat->N, 1, mat->rhs_data, mat->N, SLU_DN, SLU_D, SLU_GE);
  mat->X_data = polymec_malloc(sizeof(real_t) * mat->N);
  dCreate_Dense_Matrix(&mat->X, mat->N, 1, mat->X_data, mat->N, SLU_DN, SLU_D, SLU_GE);
  StatInit(&mat->stat);
  mat->cperm = NULL;
  mat->rperm = NULL;
  set_default_options(&mat->options);
  mat->options.ColPerm = NATURAL;
  mat->options.Fact = DOFACT;
  mat->etree = NULL;

  char name[1024];
  snprintf(name, 1024, "Sparse local matrix (N = %d)", mat->N);
  local_matrix_vtable vtable = {.dtor = slm_dtor,
                                .zero = slm_zero,
                                .add_identity = slm_add_identity,
                                .add_column_vector = slm_add_column_vector,
                                .solve = slm_solve,
                                .fprintf = slm_fprintf};
  return local_matrix_new(name, mat, vtable);
}

