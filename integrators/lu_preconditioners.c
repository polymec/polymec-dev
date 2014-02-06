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

#include <gc/gc.h>
#include "slu_ddefs.h"
#include "slu_util.h"
#include "core/array_utils.h"
#include "core/sundials_helpers.h"
#include "integrators/lu_preconditioners.h"

typedef struct 
{
  adj_graph_t* sparsity;
  adj_graph_coloring_t* coloring;
  int (*F)(void* context, real_t t, real_t* x, real_t* F);
  void (*communicate)(void* context, real_t t, real_t* x);
  void* context;

  int N; // Number of rows in matrix.

  // Work vectors.
  int num_work_vectors;
  real_t** work;

  SuperMatrix rhs, L, U;
  int *rperm, *cperm;
  superlu_options_t options;
  SuperLUStat_t stat;

  // ILU parameters (if any).
  ilu_params_t* ilu_params;
} lu_preconditioner_t;

static void supermatrix_scale_and_shift(void* context, real_t gamma)
{
  SuperMatrix* mat = context;
  int num_cols = mat->ncol;
  NCformat* data = mat->Store;
  real_t* Aij = data->nzval;
  for (int j = 0; j < num_cols; ++j)
  {
    int col_index = data->colptr[j];

    // Scale and shift diagonal value.
    Aij[col_index] = gamma * Aij[data->colptr[j]] + 1.0;

    // Scale off-diagonal values.
    size_t num_rows = data->colptr[j+1] - col_index;
    for (int i = 1; i < num_rows; ++i)
      Aij[col_index + i] *= gamma;
  }
}

static real_t supermatrix_coeff(void* context, int i, int j)
{
  SuperMatrix* mat = context;
  NCformat* data = mat->Store;
  real_t* Aij = data->nzval;
  if (i == j)
    return Aij[data->colptr[i]];
  else
  {
    int col_index = data->colptr[j];
    size_t num_rows = data->colptr[j+1] - col_index;
    int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, i);
    if (entry == NULL)
      return 0.0;
    else
    {
      size_t offset = entry - &data->rowind[col_index];
      return Aij[data->colptr[j] + offset];
    }
  }
}

static void supermatrix_fprintf(void* context, FILE* stream)
{
  SuperMatrix* A = context;
  fprintf(stream, "\nCompCol matrix P:\n");
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


static void supermatrix_dtor(void* context)
{
  SuperMatrix* matrix = context;
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
  free(matrix);
}

static preconditioner_matrix_t* lu_preconditioner_matrix(void* context)
{
  lu_preconditioner_t* precond = context;

  SuperMatrix* A = malloc(sizeof(SuperMatrix));
  
  // Fetch sparsity information from the graph.
  adj_graph_t* graph = precond->sparsity;
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

  preconditioner_matrix_vtable vtable = {.scale_and_shift = supermatrix_scale_and_shift,
                                         .coeff = supermatrix_coeff,
                                         .fprintf = supermatrix_fprintf,
                                         .dtor = supermatrix_dtor};
  return preconditioner_matrix_new("SuperMatrix", A, vtable, num_rows);
}

// Here's our finite difference implementation of the Jacobian matrix-vector 
// product. 
static void finite_diff_Jv(int (*F)(void* context, real_t t, real_t* x, real_t* F), 
                           void (*communicate)(void* context, real_t t, real_t* x),
                           void* context, 
                           real_t* x, 
                           real_t t, 
                           int num_rows,
                           real_t* v, 
                           real_t** work, 
                           real_t* Jv)
{
  real_t eps = rsqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(x).
  // work[2] == u + eps*v
  // work[3] == F(x + eps*v)

  // u + eps*v -> work[2].
  for (int i = 0; i < num_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v) -> work[3].
  if (communicate != NULL)
    communicate(context, t, work[2]);
  F(context, t, work[2], work[3]);

  // (F(x + eps*v) - F(x)) / eps -> Jv
  for (int i = 0; i < num_rows; ++i)
    Jv[i] = (work[3][i] - work[1][i]) / eps;
}

static void insert_Jv_into_supermatrix(adj_graph_t* graph, 
                                       adj_graph_coloring_t* coloring, 
                                       int color, 
                                       real_t* Jv, 
                                       SuperMatrix* J)
{
  NCformat* data = J->Store;
  real_t* Jij = data->nzval;
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    if (Jv[i] != 0.0) // <--- FIXME: Is this correct??
    {
      // Fill in the diagonal element.
      Jij[data->colptr[i]] = Jv[i];
      
      // Fill in off-diagonal column values.
      int pos = 0, j;
      while (adj_graph_next_edge(graph, i, &pos, &j))
      {
        int col_index = data->colptr[i];
        size_t num_rows = data->colptr[i+1] - col_index;
        int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, j);
        ASSERT(entry != NULL);
        size_t offset = entry - &data->rowind[col_index];
        Jij[data->colptr[i] + offset] = Jv[j];
      }
    }
  }
}

static void lu_preconditioner_compute_jacobian(void* context, real_t t, real_t* x, preconditioner_matrix_t* mat)
{
  lu_preconditioner_t* precond = context;
  adj_graph_t* graph = precond->sparsity;
  adj_graph_coloring_t* coloring = precond->coloring;
  real_t** work = precond->work;

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int num_rows = adj_graph_num_vertices(graph);
  real_t* Jv = malloc(sizeof(real_t) * num_rows);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * num_rows);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      work[0][i] = 1.0;

    // We evaluate F(x) and place it into work[1].
    if (precond->communicate != NULL)
      precond->communicate(precond->context, t, x);
    precond->F(precond->context, t, x, work[1]);

    // Now evaluate the matrix-vector product.
    memset(Jv, 0, sizeof(real_t) * num_rows);
    finite_diff_Jv(precond->F, precond->communicate, precond->context, x, t, num_rows, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    SuperMatrix* J = preconditioner_matrix_context(mat);
    insert_Jv_into_supermatrix(graph, coloring, c, Jv, J);
  }
  free(Jv);
}

static bool lu_preconditioner_solve(void* context, preconditioner_matrix_t* A, real_t* B)
{
  lu_preconditioner_t* precond = context;
  SuperMatrix* mat = preconditioner_matrix_context(A);

  // Copy B to the rhs vector.
  DNformat* rhs = precond->rhs.Store;
  memcpy(rhs->nzval, B, sizeof(real_t) * precond->N);

  if (precond->cperm == NULL)
  {
    precond->cperm = intMalloc(precond->N);
    precond->rperm = intMalloc(precond->N);
  }

  // Do the solve.
  int info;
  dgssv(&precond->options, mat, precond->cperm, precond->rperm, 
        &precond->L, &precond->U, &precond->rhs, &precond->stat, &info);
  // SuperLU claims that it can re-use factorization info with the same 
  // sparsity pattern, but this appears not to be the case at the moment.
  // precond->options.Fact = SamePattern;

  bool success = (info == 0);

  // Copy the rhs vector to B.
  memcpy(B, rhs->nzval, sizeof(real_t) * precond->N);

  return success;
}

static void lu_preconditioner_dtor(void* context)
{
  lu_preconditioner_t* precond = context;
  if (precond->cperm != NULL)
  {
    SUPERLU_FREE(precond->cperm);
    SUPERLU_FREE(precond->rperm);
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
    Destroy_SuperMatrix_Store(&precond->rhs);
  }
  for (int i = 0; i < precond->num_work_vectors; ++i)
    free(precond->work[i]);
  free(precond->work);
  adj_graph_coloring_free(precond->coloring);
  precond->ilu_params = NULL;
  free(precond);
}

preconditioner_t* lu_preconditioner_new(void* context,
                                        int (*residual_func)(void* context, real_t t, real_t* x, real_t* F),
                                        void (*communication_func)(void* context, real_t t, real_t* x),
                                        adj_graph_t* sparsity)
{
  lu_preconditioner_t* precond = malloc(sizeof(lu_preconditioner_t));
  precond->sparsity = sparsity;
  precond->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
  precond->F = residual_func;
  precond->communicate = communication_func;
  precond->context = context;
  precond->ilu_params = NULL;

  // Preconditioner data.
  precond->N = adj_graph_num_vertices(precond->sparsity);
  real_t* rhs = malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->rhs, precond->N, 1, rhs, precond->N, SLU_DN, SLU_D, SLU_GE);
  StatInit(&precond->stat);
  precond->cperm = NULL;
  precond->rperm = NULL;
  set_default_options(&precond->options);
  precond->options.ColPerm = NATURAL;
  precond->options.Fact = DOFACT;

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = malloc(sizeof(real_t) * precond->N);

  preconditioner_vtable vtable = {.matrix = lu_preconditioner_matrix,
                                  .compute_jacobian = lu_preconditioner_compute_jacobian,
                                  .solve = lu_preconditioner_solve,
                                  .dtor = lu_preconditioner_dtor};
  return preconditioner_new("LU preconditioner", precond, vtable);
}
                                        
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
  params->variant = ILU_SILU;
  params->fill_tolerance = 0.01;
  return params;
}

static bool ilu_preconditioner_solve(void* context, preconditioner_matrix_t* A, real_t* B)
{
  // FIXME
  return true;
}

// ILU preconditioner.
preconditioner_t* ilu_preconditioner_new(void* context,
                                         int (*residual_func)(void* context, real_t t, real_t* x, real_t* F),
                                         void (*communication_func)(void* context, real_t t, real_t* x),
                                         adj_graph_t* sparsity, 
                                         ilu_params_t* ilu_params)
{
  lu_preconditioner_t* precond = malloc(sizeof(lu_preconditioner_t));
  precond->sparsity = sparsity;
  precond->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
  precond->F = residual_func;
  precond->communicate = communication_func;
  precond->context = context;
  precond->ilu_params = ilu_params;

  // Make work vectors.
  precond->N = adj_graph_num_vertices(precond->sparsity);
  precond->num_work_vectors = 4;
  precond->work = malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = malloc(sizeof(real_t) * precond->N);

  preconditioner_vtable vtable = {.matrix = lu_preconditioner_matrix,
                                  .compute_jacobian = lu_preconditioner_compute_jacobian,
                                  .solve = ilu_preconditioner_solve,
                                  .dtor = lu_preconditioner_dtor};
  return preconditioner_new("ILU preconditioner", precond, vtable);
}

