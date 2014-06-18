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
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*dae_F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* context;

  int N; // Number of rows in matrix.

  // Work vectors.
  int num_work_vectors;
  real_t** work;
  int* etree;

  SuperMatrix rhs, X, L, U;
  real_t *rhs_data, *X_data;
  int *rperm, *cperm;
  superlu_options_t options;
  SuperLUStat_t stat;

  // Preconditioner matrix.
  SuperMatrix* P;

  // ILU parameters (if any).
  ilu_params_t* ilu_params;
} lu_preconditioner_t;

static void supermatrix_eye(SuperMatrix* mat, real_t diag_val)
{
  int num_cols = mat->ncol;
  NCformat* data = mat->Store;
  real_t* Aij = data->nzval;

  // Zero the matrix coefficients.
  memset(Aij, 0, data->colptr[num_cols]);

  // Set the diagonal values.
  for (int j = 0; j < num_cols; ++j)
  {
    int col_index = data->colptr[j];
    Aij[col_index] = diag_val;
  }
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

// Here's our finite difference implementation of the dF/dx matrix-vector 
// product. 
static void finite_diff_dFdx_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval), 
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
  // work[3] == F(t, x + eps*v, xdot)

  // x + eps*v -> work[2].
  for (int i = 0; i < num_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v, xdot) -> work[3].
  F(context, t, work[2], xdot, work[3]);

  // (F(t, x + eps*v, xdot) - F(t, x, xdot)) / eps -> dFdx_v
  for (int i = 0; i < num_rows; ++i)
    dFdx_v[i] = (work[3][i] - work[1][i]) / eps;
}

// Here's the same matrix-vector product for dF/d(xdot).  
static void finite_diff_dFdxdot_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval), 
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
    work[2][i] = x[i] + eps*v[i];

  // F(t, x, xdot + eps*v) -> work[3].
  F(context, t, x, work[2], work[3]);

  // (F(t, x, xdot + eps*v, xdot) - F(t, x, xdot)) / eps -> dFdxdot_v
  for (int i = 0; i < num_rows; ++i)
    dFdxdot_v[i] = (work[3][i] - work[1][i]) / eps;
}

static void add_Jv_into_supermatrix(adj_graph_t* graph, 
                                    adj_graph_coloring_t* coloring, 
                                    int color, 
                                    real_t factor, 
                                    real_t* Jv, 
                                    SuperMatrix* J)
{
  NCformat* data = J->Store;
  real_t* Jij = data->nzval;
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    // Add in the diagonal element.
    Jij[data->colptr[i]] += Jv[i];
      
    // Add in off-diagonal column values.
    int pos = 0, j;
    while (adj_graph_next_edge(graph, i, &pos, &j))
    {
      int col_index = data->colptr[i];
      size_t num_rows = data->colptr[i+1] - col_index;
      int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, j);
      ASSERT(entry != NULL);
      size_t offset = entry - &data->rowind[col_index];
      Jij[data->colptr[i] + offset] += factor * Jv[j];
    }
  }
}

// This function adapts non-DAE functions F(t, x) to DAE ones F(t, x, xdot).
static int F_adaptor(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval)
{
  ASSERT(xdot == NULL);

  // We are passed the actual preconditioner as our context pointer, so get the 
  // "real" one here.
  lu_preconditioner_t* pc = context;
  void* ctx = pc->context;
  return pc->F(ctx, t, x, Fval);
}

static void lu_preconditioner_setup(void* context, real_t alpha, real_t beta, real_t gamma,
                                    real_t t, real_t* x, real_t* xdot)
{
  lu_preconditioner_t* precond = context;
  adj_graph_t* graph = precond->sparsity;
  adj_graph_coloring_t* coloring = precond->coloring;
  real_t** work = precond->work;

  // Zero our preconditioner matrix and set the diagonal term.
  supermatrix_eye(precond->P, alpha);

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
  real_t* Jv = polymec_malloc(sizeof(real_t) * num_rows);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * num_rows);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      work[0][i] = 1.0;

    // We evaluate F(x) and place it into work[1].
    F(F_context, t, x, xdot, work[1]);

    // Now evaluate the matrix-vector product for dF/dx and stash it.
    memset(Jv, 0, sizeof(real_t) * num_rows);
    finite_diff_dFdx_v(F, F_context, t, x, xdot, num_rows, work[0], work, Jv);
    add_Jv_into_supermatrix(graph, coloring, c, beta, Jv, precond->P);

    // Do the same for dF/d(xdot).
    if ((gamma != 0.0) && (xdot != NULL))
    {
      memset(Jv, 0, sizeof(real_t) * num_rows);
      finite_diff_dFdxdot_v(F, F_context, t, x, xdot, num_rows, work[0], work, Jv);
      add_Jv_into_supermatrix(graph, coloring, c, gamma, Jv, precond->P);
    }
  }
  polymec_free(Jv);
}

static bool lu_preconditioner_solve(void* context, real_t* B)
{
  lu_preconditioner_t* precond = context;
  SuperMatrix* mat = precond->P;

  // Copy B to the rhs vector.
  DNformat* rhs = precond->rhs.Store;
  memcpy(rhs->nzval, B, sizeof(real_t) * precond->N);

  if (precond->cperm == NULL)
  {
    precond->cperm = intMalloc(precond->N);
    precond->rperm = intMalloc(precond->N);
  }
  else if (precond->options.Fact != SamePattern)
  {
    // Jettison the existing factorization.
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
  }

  // Do the solve.
  int info;
  dgssv(&precond->options, mat, precond->cperm, precond->rperm, 
        &precond->L, &precond->U, &precond->rhs, &precond->stat, &info);
  // SuperLU claims that it can re-use factorization info with the same 
  // sparsity pattern, but this appears not to be the case at the moment.
  // precond->options.Fact = SamePattern;

  bool success = (info == 0);

  if (success)
  {
    // Copy the rhs vector to B.
    memcpy(B, rhs->nzval, sizeof(real_t) * precond->N);
  }

  return success;
}

static void lu_preconditioner_fprintf(void* context, FILE* stream)
{
  lu_preconditioner_t* precond = context;
  SuperMatrix* A = precond->P;
  fprintf(stream, "\nCompCol preconditioner matrix P:\n");
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

static void lu_preconditioner_dtor(void* context)
{
  lu_preconditioner_t* precond = context;
  if (precond->cperm != NULL)
  {
    SUPERLU_FREE(precond->cperm);
    SUPERLU_FREE(precond->rperm);
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
  }
  Destroy_SuperMatrix_Store(&precond->rhs);
  polymec_free(precond->rhs_data);
  Destroy_SuperMatrix_Store(&precond->X);
  supermatrix_free(precond->P);
  polymec_free(precond->X_data);
  StatFree(&precond->stat);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    polymec_free(precond->work[i]);
  polymec_free(precond->work);
  if (precond->etree != NULL)
    polymec_free(precond->etree);
  adj_graph_coloring_free(precond->coloring);
  precond->ilu_params = NULL;
  polymec_free(precond);
}

static preconditioner_t* general_lu_preconditioner_new(void* context,
                                                       adj_graph_t* sparsity)
{
  lu_preconditioner_t* precond = polymec_malloc(sizeof(lu_preconditioner_t));
  precond->sparsity = sparsity;
  precond->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
  log_debug("LU preconditioner: graph coloring produced %d colors.", 
            adj_graph_coloring_num_colors(precond->coloring));
  precond->context = context;
  precond->ilu_params = NULL;
  precond->P = supermatrix_new(sparsity);

  // Preconditioner data.
  precond->N = adj_graph_num_vertices(precond->sparsity);
  precond->rhs_data = polymec_malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->rhs, precond->N, 1, precond->rhs_data, precond->N, SLU_DN, SLU_D, SLU_GE);
  precond->X_data = polymec_malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->X, precond->N, 1, precond->X_data, precond->N, SLU_DN, SLU_D, SLU_GE);
  StatInit(&precond->stat);
  precond->cperm = NULL;
  precond->rperm = NULL;
  set_default_options(&precond->options);
  precond->options.ColPerm = NATURAL;
  precond->options.Fact = DOFACT;
  precond->etree = NULL;
  precond->F = NULL;
  precond->dae_F = NULL;

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = polymec_malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = polymec_malloc(sizeof(real_t) * precond->N);

  preconditioner_vtable vtable = {.setup = lu_preconditioner_setup,
                                  .solve = lu_preconditioner_solve,
                                  .fprintf = lu_preconditioner_fprintf,
                                  .dtor = lu_preconditioner_dtor};
  return preconditioner_new("LU preconditioner", precond, vtable);
}
   
preconditioner_t* lu_preconditioner_new(void* context,
                                        int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                        adj_graph_t* sparsity)
{
  ASSERT(F != NULL);
  preconditioner_t* precond = general_lu_preconditioner_new(context, sparsity);
  lu_preconditioner_t* lu = preconditioner_context(precond);
  lu->F = F;
  return precond;
}

preconditioner_t* lu_dae_preconditioner_new(void* context,
                                            int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                            adj_graph_t* sparsity)
{
  ASSERT(F != NULL);
  preconditioner_t* precond = general_lu_preconditioner_new(context, sparsity);
  lu_preconditioner_t* lu = preconditioner_context(precond);
  lu->dae_F = F;
  return precond;
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
  params->milu_variant = ILU_SILU;
  params->fill_tolerance = 0.01;
  params->norm = ILU_LINF;
  return params;
}

static bool ilu_preconditioner_solve(void* context, real_t* B)
{
  lu_preconditioner_t* precond = context;
  SuperMatrix* mat = precond->P;

  // Copy B to the rhs vector.
  DNformat* rhs = precond->rhs.Store;
  memcpy(rhs->nzval, B, sizeof(real_t) * precond->N);

  if (precond->cperm == NULL)
  {
    precond->cperm = intMalloc(precond->N);
    precond->rperm = intMalloc(precond->N);
  }
  else if (precond->options.Fact != SamePattern)
  {
    // Jettison the existing factorization.
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
  }


  // Do the (approximate) solve.
  int info;
  char equed;
  int lwork = 0; // Indicates SuperLU internal memory allocation.
  real_t R[precond->N], C[precond->N], recip_pivot_growth, cond_number;
  mem_usage_t mem_usage;
  dgsisx(&precond->options, mat, precond->cperm, precond->rperm, 
         precond->etree, &equed, R, C, &precond->L, &precond->U, 
         NULL, lwork, &precond->rhs, &precond->X, &recip_pivot_growth, &cond_number,
         &mem_usage, &precond->stat, &info);

  bool success = (info == 0);

  // Copy the rhs vector to B.
  memcpy(B, rhs->nzval, sizeof(real_t) * precond->N);

  return success;
}

// ILU preconditioner.
static preconditioner_t* general_ilu_preconditioner_new(void* context,
                                                        adj_graph_t* sparsity, 
                                                        ilu_params_t* ilu_params)
{
  lu_preconditioner_t* precond = polymec_malloc(sizeof(lu_preconditioner_t));
  precond->sparsity = sparsity;
  precond->coloring = adj_graph_coloring_new(sparsity, SMALLEST_LAST);
  log_debug("ILU preconditioner: graph coloring produced %d colors.", 
            adj_graph_coloring_num_colors(precond->coloring));
  precond->context = context;
  precond->P = supermatrix_new(sparsity);

  // Copy options into place.
  ilu_set_default_options(&precond->options);
  precond->ilu_params = ilu_params;
  precond->options.ILU_DropRule = ilu_params->drop_rule;
  precond->options.ILU_DropTol = ilu_params->drop_tolerance;
  precond->options.ILU_FillFactor = ilu_params->fill_factor;
  precond->options.ILU_FillTol = ilu_params->fill_tolerance;
  if (ilu_params->norm == ILU_L1)
    precond->options.ILU_Norm = ONE_NORM;
  else if (ilu_params->norm == ILU_L2)
    precond->options.ILU_Norm = TWO_NORM;
  else
    precond->options.ILU_Norm = INF_NORM;
  if (ilu_params->milu_variant == ILU_SILU)
    precond->options.ILU_MILU = SILU;
  else if (ilu_params->milu_variant == ILU_MILU1)
    precond->options.ILU_MILU = SMILU_1;
  else if (ilu_params->milu_variant == ILU_MILU2)
    precond->options.ILU_MILU = SMILU_2;
  else 
    precond->options.ILU_MILU = SMILU_3;
  precond->options.ILU_MILU_Dim = 3;
  if (ilu_params->row_perm == ILU_NO_ROW_PERM)
    precond->options.RowPerm = NOROWPERM;
  else
    precond->options.RowPerm = LargeDiag;

  // Preconditioner data.
  precond->N = adj_graph_num_vertices(precond->sparsity);
  precond->rhs_data = polymec_malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->rhs, precond->N, 1, precond->rhs_data, precond->N, SLU_DN, SLU_D, SLU_GE);
  precond->X_data = polymec_malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->X, precond->N, 1, precond->X_data, precond->N, SLU_DN, SLU_D, SLU_GE);
  StatInit(&precond->stat);
  precond->cperm = NULL;
  precond->rperm = NULL;
  precond->options.ColPerm = NATURAL;
  precond->options.Fact = DOFACT;
  precond->etree = polymec_malloc(sizeof(int) * precond->N);
  precond->F = NULL;
  precond->dae_F = NULL;

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = polymec_malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = polymec_malloc(sizeof(real_t) * precond->N);

  preconditioner_vtable vtable = {.setup = lu_preconditioner_setup,
                                  .solve = ilu_preconditioner_solve,
                                  .fprintf = lu_preconditioner_fprintf,
                                  .dtor = lu_preconditioner_dtor};
  return preconditioner_new("ILU preconditioner", precond, vtable);
}

preconditioner_t* ilu_preconditioner_new(void* context,
                                         int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                         adj_graph_t* sparsity,
                                         ilu_params_t* ilu_params)
{
  ASSERT(F != NULL);
  preconditioner_t* precond = general_ilu_preconditioner_new(context, sparsity, ilu_params);
  lu_preconditioner_t* ilu = preconditioner_context(precond);
  ilu->F = F;
  return precond;
}

preconditioner_t* ilu_dae_preconditioner_new(void* context,
                                             int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                             adj_graph_t* sparsity,
                                             ilu_params_t* ilu_params)
{
  ASSERT(F != NULL);
  preconditioner_t* precond = general_ilu_preconditioner_new(context, sparsity, ilu_params);
  lu_preconditioner_t* ilu = preconditioner_context(precond);
  ilu->dae_F = F;
  return precond;
}


