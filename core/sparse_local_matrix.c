// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "slu_ddefs.h"
#include "slu_util.h"
#include "core/array_utils.h"
#include "core/sparse_local_matrix.h"

typedef struct 
{
  int N; // Number of rows in matrix.

  // Personal copy of the sparsity graph.
  adj_graph_t* sparsity;

  // Matrix data.
  int* etree;
  char equil;
  SuperMatrix rhs, X, L, U;
  double *rhs_data, *X_data;
  double *R, *C;
  int *rperm, *cperm;
  superlu_options_t options;
  SuperLUStat_t stat;
  SuperMatrix* A; // <-- matrix itself!.

  // ILU parameters (if any).
  ilu_params_t* ilu_params;
} slm_t;

// Globals borrowed from SuperLU.
const int ILU_DROP_BASIC = DROP_BASIC;
const int ILU_DROP_PROWS = DROP_PROWS;
const int ILU_DROP_COLUMN = DROP_COLUMN;
const int ILU_DROP_AREA = DROP_AREA;
const int ILU_DROP_DYNAMIC = DROP_DYNAMIC;
const int ILU_DROP_INTERP = DROP_INTERP;

ilu_params_t* ilu_params_new()
{
  ilu_params_t* params = GC_MALLOC(sizeof(ilu_params_t));
  params->diag_pivot_threshold = 0.1;
  params->row_perm = ILU_LARGE_DIAG_PERM;
  params->drop_rule = ILU_DROP_BASIC | ILU_DROP_AREA;
  params->drop_tolerance = 1e-4;
  params->fill_factor = 10.0;
  params->milu_variant = ILU_SILU;
  params->fill_tolerance = 0.01;
  params->norm = ILU_LINF;
  return params;
}

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

static void* slm_clone(void* context)
{
  slm_t* mat = context;
  slm_t* clone = polymec_malloc(sizeof(slm_t));
  clone->sparsity = adj_graph_clone(mat->sparsity);
  clone->A = supermatrix_new(clone->sparsity);
  int nnz = ((NCformat*)clone->A->Store)->nnz;
  real_t* Aij = ((NCformat*)clone->A->Store)->nzval;
  real_t* Bij = ((NCformat*)mat->A->Store)->nzval;
  memcpy(Aij, Bij, sizeof(real_t) * nnz);

  clone->N = mat->N;
  clone->rhs_data = polymec_malloc(sizeof(double) * clone->N);
  memcpy(clone->rhs_data, mat->rhs_data, sizeof(double) * clone->N);
  dCreate_Dense_Matrix(&clone->rhs, clone->N, 1, clone->rhs_data, clone->N, SLU_DN, SLU_D, SLU_GE);
  clone->X_data = polymec_malloc(sizeof(double) * clone->N);
  memcpy(clone->X_data, mat->X_data, sizeof(double) * clone->N);
  dCreate_Dense_Matrix(&clone->X, clone->N, 1, clone->X_data, clone->N, SLU_DN, SLU_D, SLU_GE);
  clone->R = polymec_malloc(sizeof(double) * clone->N);
  memcpy(clone->R, mat->R, sizeof(double) * clone->N);
  clone->C = polymec_malloc(sizeof(double) * clone->N);
  memcpy(clone->C, mat->C, sizeof(double) * clone->N);
  StatInit(&clone->stat);

  clone->cperm = NULL;
  clone->rperm = NULL;
  clone->options = mat->options;
  clone->options.Fact = DOFACT;
  clone->etree = NULL;

  return clone;
}

static void slm_zero(void* context)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;
  NCformat* data = A->Store;
  real_t* Aij = data->nzval;
  memset(Aij, 0, sizeof(real_t) * data->nnz);

  // We have to redo the factorization.
  mat->options.Fact = DOFACT;
}

static int slm_num_columns(void* context, int row)
{
  slm_t* mat = context;
  return 1 + adj_graph_num_edges(mat->sparsity, row);
}

static void slm_get_columns(void* context, int row, int* columns)
{
  slm_t* mat = context;
  columns[0] = row;
  memcpy(&columns[1], adj_graph_edges(mat->sparsity, row), 
         sizeof(int) * adj_graph_num_edges(mat->sparsity, row));
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

  // We have to redo the factorization.
  mat->options.Fact = DOFACT;
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

  if (column >= mat->N) return;
  int col_index = data->colptr[column];

  // Add in the row values.
  size_t num_rows = data->colptr[column+1] - col_index;
  for (int r = 0; r < num_rows; ++r)
  {
    int row = data->rowind[col_index+r];
    Jij[col_index+r] += scale_factor * column_vector[row];
  }

  // We have to redo the factorization.
  mat->options.Fact = DOFACT;
}

static void slm_add_row_vector(void* context, 
                               real_t scale_factor,
                               int row,
                               real_t* row_vector)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;

  NCformat* data = A->Store;
  real_t* Aij = data->nzval;

  int i = row;

  if (i >= mat->N) return;

  // Add in the diagonal element.
  Aij[data->colptr[i]] += scale_factor * row_vector[i];
      
  // Add in off-diagonal column values. 
  int pos = 0, j;
  while (adj_graph_next_edge(mat->sparsity, i, &pos, &j))
  {
    int col_index = data->colptr[j];
    size_t num_rows = data->colptr[j+1] - col_index;
    if (num_rows > 1)
    {
      int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, i);
      ASSERT(entry != NULL);
      size_t offset = entry - &data->rowind[col_index];
      Aij[data->colptr[j] + offset] += scale_factor * row_vector[i];
    }
  }

  // We have to redo the factorization.
  mat->options.Fact = DOFACT;
}

static bool slm_solve(void* context, real_t* B, real_t* x)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;

  // Copy B to the rhs vector.
  DNformat* rhs = mat->rhs.Store;
  double* Bi = rhs->nzval;
  for (int i = 0; i < mat->N; ++i)
    Bi[i] = (double)B[i];

  if (mat->cperm == NULL)
  {
    mat->cperm = intMalloc(mat->N);
    mat->rperm = intMalloc(mat->N);
  }
  else if (mat->options.Fact != FACTORED)
  {
    // Jettison the existing factorization.
    Destroy_SuperNode_Matrix(&mat->L);
    Destroy_CompCol_Matrix(&mat->U);
  }

  // Do the solve.
  int info = 0, lwork = 0;
  void* work = NULL;
  double ferr, berr;
  GlobalLU_t glu; // "Global data structure" for SuperLU for helping with factorizations.
  double recip_pivot_growth = 1.0, rcond = 1.0;
  mem_usage_t mem_usage;
  if (mat->ilu_params == NULL)
  {
    // Factorize if necessary.
    if (mat->options.Fact == DOFACT)
    {
      int rhs_ncol = mat->rhs.ncol;
      mat->rhs.ncol = 0;
      polymec_suspend_fpe();
      dgssvx(&mat->options, A, mat->cperm, mat->rperm, mat->etree, &mat->equil, 
             mat->R, mat->C, &mat->L, &mat->U, work, lwork, &mat->rhs, &mat->X, 
             &recip_pivot_growth, &rcond, &ferr, &berr, &glu, &mem_usage, &mat->stat, &info);
      polymec_restore_fpe();
      mat->rhs.ncol = rhs_ncol;

      if ((info == 0) || (info == A->ncol+1))
      {
        if (mat->equil != 'N')
        {
          if (mat->equil == 'R')
            log_debug("slm_solve: performed row equilibration.");
          else if (mat->equil == 'C')
            log_debug("slm_solve: performed column equilibration.");
          else if (mat->equil == 'B')
            log_debug("slm_solve: performed row/column equilibration.");
        }
        log_debug("slm_solve: L has %d nonzeros, U has %d.", 
                  ((SCformat*)mat->L.Store)->nnz, ((NCformat*)mat->U.Store)->nnz);
#ifndef NDEBUG
        log_debug("slm_solve: recip pivot growth = %g, condition number = %g.", 
                  recip_pivot_growth, rcond);
        if (recip_pivot_growth < 0.01)
        {
          log_detail("slm_solve: WARNING: recip pivot growth for ILU factorization << 1.");
          log_detail("slm_solve: WARNING: Stability of LU factorization could be poor.");
        }
#endif

        // Reuse the factorization.
        mat->options.Fact = FACTORED;
      }
      else
        log_debug("slm_solve: LU factorization failed.");
    }

    // Solve the factored system.
    if ((info == 0) || (info == A->ncol+1))
    {
      polymec_suspend_fpe();
      dgssvx(&mat->options, A, mat->cperm, mat->rperm, mat->etree, &mat->equil, 
             mat->R, mat->C, &mat->L, &mat->U, work, lwork, &mat->rhs, 
             &mat->X, &recip_pivot_growth, &rcond, &ferr, &berr, &glu, &mem_usage, 
             &mat->stat, &info);
      polymec_restore_fpe();
    }
  }
  else
  {
    // Incomplete LU factorization.

    // Factorize if necessary.
    if (mat->options.Fact == DOFACT)
    {
      int rhs_ncol = mat->rhs.ncol;
      mat->rhs.ncol = 0;
      polymec_suspend_fpe();
      dgsisx(&mat->options, A, mat->cperm, mat->rperm, mat->etree, &mat->equil, 
             mat->R, mat->C, &mat->L, &mat->U, NULL, 0, &mat->rhs, &mat->X, 
             &recip_pivot_growth, &rcond, &glu, &mem_usage, &mat->stat, &info);
      polymec_restore_fpe();
      mat->rhs.ncol = rhs_ncol;

      if ((info == 0) || (info == A->ncol+1))
      {
        if (mat->equil != 'N')
        {
          if (mat->equil == 'R')
            log_debug("slm_solve: performed row equilibration.");
          else if (mat->equil == 'C')
            log_debug("slm_solve: performed column equilibration.");
          else if (mat->equil == 'B')
            log_debug("slm_solve: performed row/column equilibration.");
        }
#ifndef NDEBUG
        log_debug("slm_solve: recip pivot growth = %g, condition number = %g.", 
                  recip_pivot_growth, rcond);
        if (recip_pivot_growth < 0.01)
        {
          log_detail("slm_solve: WARNING: recip pivot growth for ILU factorization << 1.");
          log_detail("slm_solve: WARNING: Stability of LU factorization could be poor.");
        }
#endif

        // Reuse the factorization.
        mat->options.Fact = FACTORED;
      }
      else
        log_debug("slm_solve: incomplete LU factorization failed.");
    }

    // Solve the factored system.
    if ((info == 0) || (info == A->ncol+1))
    {
      polymec_suspend_fpe();
      dgsisx(&mat->options, A, mat->cperm, mat->rperm, mat->etree, &mat->equil, 
             mat->R, mat->C, &mat->L, &mat->U, NULL, 0, &mat->rhs, &mat->X, 
             &recip_pivot_growth, &rcond, &glu, &mem_usage, &mat->stat, &info);
      polymec_restore_fpe();
    }
  }

  bool success = ((info == 0) || (info == A->ncol+1));
  if (success)
  {
    // Copy the output vector to x.
    double* X = ((DNformat*)mat->X.Store)->nzval;
    for (int i = 0; i < mat->N; ++i)
      x[i] = (real_t)X[i];
  }
  else
  {
    ASSERT(info > 0);
    if (mat->ilu_params == NULL)
    {
      log_debug("slm_solve: LU solve failed.");
      log_debug("slm_solve: (U is singular: U(%d, %d) = 0.)", info-1, info-1);
    }
    else 
    {
      log_debug("slm_solve: ILU solve failed.");
      if (info < A->ncol)
        log_debug("slm_solve: (number of zero pivots in U = %d.)", info);
      else if (info == (A->ncol + 1))
        log_debug("slm_solve: (U is nonsingular but rcond = %g.)", rcond);
      else
        log_debug("slm_solve: (Memory allocation failure.)");
    }
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

static real_t slm_value(void* context, int i, int j)
{
  slm_t* mat = context;
  ASSERT(i < mat->N);
  ASSERT(j < mat->N);

  SuperMatrix* A = mat->A;

  NCformat* data = A->Store;
  real_t* Aij = data->nzval;

  if (j == i) // diagonal value
    return Aij[data->colptr[i]];

  int col_index = data->colptr[j];
  size_t num_rows = data->colptr[j+1] - col_index;
  int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, i);
  if (entry == NULL)
    return 0.0;
  size_t offset = entry - &data->rowind[col_index];
  return Aij[data->colptr[j] + offset];
}

static void slm_set_value(void* context, int i, int j, real_t value)
{
  slm_t* mat = context;
  ASSERT(i < mat->N);
  ASSERT(j < mat->N);

  SuperMatrix* A = mat->A;

  NCformat* data = A->Store;
  real_t* Aij = data->nzval;

  if (j == i) // diagonal value
    Aij[data->colptr[i]] = value;
  else
  {
    int col_index = data->colptr[j];
    size_t num_rows = data->colptr[j+1] - col_index;
    int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, i);
    if (entry != NULL)
    {
      size_t offset = entry - &data->rowind[col_index];
      Aij[data->colptr[j] + offset] = value;
    }
  }

  // We have to redo the factorization.
  mat->options.Fact = DOFACT;
}

static void slm_get_diag(void* context, real_t* diag)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;
  NCformat* data = A->Store;
  real_t* Aij = data->nzval;

  for (int i = 0; i < mat->N; ++i)
    diag[i] = Aij[data->colptr[i]];
}

static void slm_matvec(void* context, real_t* x, real_t* Ax)
{
  slm_t* mat = context;
  SuperMatrix* A = mat->A;

  NCformat* data = A->Store;
  real_t* Aij = data->nzval;

  for (int i = 0; i < mat->N; ++i)
  {
    // Add in the diagonal element.
    Ax[i] = Aij[data->colptr[i]] * x[i];
      
    // Add in off-diagonal column values. 
    int pos = 0, j;
    while (adj_graph_next_edge(mat->sparsity, i, &pos, &j))
    {
      int col_index = data->colptr[j];
      size_t num_rows = data->colptr[j+1] - col_index;
      if (num_rows > 1)
      {
        int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, i);
        ASSERT(entry != NULL);
        size_t offset = entry - &data->rowind[col_index];
        Ax[i] += Aij[data->colptr[j] + offset] * x[j];
      }
    }
  }
}

static void slm_add_matrix(void* context, real_t scale_factor, void* B)
{
  slm_t* Amat = context;
  int nnz = ((NCformat*)Amat->A->Store)->nnz;
  real_t* Aij = ((NCformat*)Amat->A->Store)->nzval;
  slm_t* Bmat = context;
  ASSERT(nnz == ((NCformat*)Bmat->A->Store)->nnz);
  real_t* Bij = ((NCformat*)Bmat->A->Store)->nzval;
  for (int i = 0; i < nnz; ++i)
    Aij[i] += scale_factor * Bij[i];

  // We have to redo the factorization.
  Amat->options.Fact = DOFACT;
}

extern double dlangs(char* norm, SuperMatrix* A); // SuperLU norm function.

static real_t slm_norm(void* context, char n)
{
  slm_t* mat = context;
  return (real_t)dlangs(&n, mat->A);
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
  supermatrix_free(mat->A);
  Destroy_SuperMatrix_Store(&mat->rhs);
  polymec_free(mat->rhs_data);
  Destroy_SuperMatrix_Store(&mat->X);
  polymec_free(mat->X_data);
  polymec_free(mat->R);
  polymec_free(mat->C);
  StatFree(&mat->stat);
  if (mat->etree != NULL)
    polymec_free(mat->etree);
  mat->ilu_params = NULL;
  adj_graph_free(mat->sparsity);
  polymec_free(mat);
}

local_matrix_t* sparse_local_matrix_new(adj_graph_t* sparsity)
{
  slm_t* mat = polymec_malloc(sizeof(slm_t));
  mat->sparsity = adj_graph_clone(sparsity); // MINE!
  mat->ilu_params = NULL;
  mat->A = supermatrix_new(sparsity);

  // Solver data.
  mat->N = adj_graph_num_vertices(sparsity);
  mat->rhs_data = polymec_malloc(sizeof(double) * mat->N);
  dCreate_Dense_Matrix(&mat->rhs, mat->N, 1, mat->rhs_data, mat->N, SLU_DN, SLU_D, SLU_GE);
  mat->X_data = polymec_malloc(sizeof(double) * mat->N);
  dCreate_Dense_Matrix(&mat->X, mat->N, 1, mat->X_data, mat->N, SLU_DN, SLU_D, SLU_GE);
  mat->R = polymec_malloc(sizeof(double) * mat->N);
  mat->C = polymec_malloc(sizeof(double) * mat->N);
  StatInit(&mat->stat);
  mat->cperm = NULL;
  mat->rperm = NULL;
  mat->etree = polymec_malloc(sizeof(int) * mat->N);
  set_default_options(&mat->options);
  mat->options.Equil = NO;
//  mat->options.ColPerm = NATURAL;
  mat->options.Fact = DOFACT;
#ifndef NDEBUG
  mat->options.PivotGrowth = YES;
  mat->options.ConditionNumber = YES;
#endif

  char name[1024];
  snprintf(name, 1024, "Sparse local matrix (N = %d)", mat->N);
  local_matrix_vtable vtable = {.clone = slm_clone,
                                .dtor = slm_dtor,
                                .zero = slm_zero,
                                .num_columns = slm_num_columns,
                                .get_columns = slm_get_columns,
                                .add_identity = slm_add_identity,
                                .add_column_vector = slm_add_column_vector,
                                .add_row_vector = slm_add_row_vector,
                                .solve = slm_solve,
                                .fprintf = slm_fprintf,
                                .value = slm_value,
                                .set_value = slm_set_value,
                                .get_diag = slm_get_diag,
                                .matvec = slm_matvec,
                                .add_matrix = slm_add_matrix,
                                .norm = slm_norm};
  return local_matrix_new(name, mat, vtable, mat->N);
}

void sparse_local_matrix_use_ilu(local_matrix_t* matrix,
                                 ilu_params_t* ilu_params)
{
  ASSERT(ilu_params != NULL);

  slm_t* mat = local_matrix_context(matrix);
  mat->ilu_params = ilu_params;
  ilu_set_default_options(&mat->options);
  mat->options.DiagPivotThresh = ilu_params->diag_pivot_threshold;
  if (ilu_params->row_perm == ILU_NO_ROW_PERM)
    mat->options.RowPerm = NOROWPERM;
  else
    mat->options.RowPerm = LargeDiag;
  mat->options.ILU_DropRule = ilu_params->drop_rule;
  mat->options.ILU_DropTol = ilu_params->drop_tolerance;
  mat->options.ILU_FillFactor = ilu_params->fill_factor;
  if (ilu_params->milu_variant == ILU_SILU)
    mat->options.ILU_MILU = SILU;
  else if (ilu_params->milu_variant == ILU_MILU1)
    mat->options.ILU_MILU = SMILU_1;
  else if (ilu_params->milu_variant == ILU_MILU2)
    mat->options.ILU_MILU = SMILU_2;
  else 
    mat->options.ILU_MILU = SMILU_3;
  mat->options.ILU_FillTol = ilu_params->fill_tolerance;
  if (ilu_params->norm == ILU_L1)
    mat->options.ILU_Norm = ONE_NORM;
  else if (ilu_params->norm == ILU_L2)
    mat->options.ILU_Norm = TWO_NORM;
  else
    mat->options.ILU_Norm = INF_NORM;
}

void sparse_local_matrix_use_lu(local_matrix_t* matrix)
{
  slm_t* mat = local_matrix_context(matrix);
  mat->ilu_params = NULL;
  set_default_options(&mat->options);
}

