// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
  SuperMatrix rhs, X, L, U;
  real_t *rhs_data, *X_data;
  int *rperm, *cperm;
  superlu_options_t options;
  SuperLUStat_t stat;
  SuperMatrix* A; // <-- matrix itself!.

  // ILU parameters (if any).
  ilu_params_t* ilu_params;
} slm_t;

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

static bool slm_solve(void* context, real_t* B)
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

  int col_index = data->colptr[i];
  size_t num_rows = data->colptr[i+1] - col_index;
  int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, j);
  if (entry == NULL)
    return 0.0;
  size_t offset = entry - &data->rowind[col_index];
  return Aij[data->colptr[i] + offset];
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

local_matrix_t* sparse_local_matrix_new(adj_graph_t* sparsity)
{
  return ilu_sparse_local_matrix_new(sparsity, NULL);
}

local_matrix_t* ilu_sparse_local_matrix_new(adj_graph_t* sparsity,
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
                                .fprintf = slm_fprintf,
                                .value = slm_value};
  return local_matrix_new(name, mat, vtable);
}

