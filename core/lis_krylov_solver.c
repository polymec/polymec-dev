// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/krylov_solver.h"

// Set up stuff for LIS.
#if POLYMEC_HAVE_MPI
#define USE_MPI 1 // needed by lis.h to enable MPI
#endif
#ifdef _LONG__DOUBLE
#undef _LONG__DOUBLE  // No long doubles, thank you!
#endif
// NOTE: We want LIS_INT == long long, which is the same size as index_t.
#ifndef _LONG__LONG
#define _LONG__LONG 
#endif
#include "lis.h"

#include "core/options.h"
#include "core/timer.h"
#include "core/array_utils.h"

// For some reason, these prototypes didn't make it into lis.h.
extern LIS_INT lis_matrix_shift_diagonal(LIS_MATRIX A, LIS_SCALAR alpha);
extern LIS_INT lis_matrix_set_destroyflag(LIS_MATRIX A, LIS_INT flag);

//------------------------------------------------------------------------
// This file implements the LIS Krylov solver.
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;
} lis_factory_t;

typedef struct
{
  LIS_SOLVER solver;
  LIS_MATRIX op;
} lis_solver_t;

typedef struct
{
  MPI_Comm comm;
  LIS_MATRIX A;
  LIS_INT block_size, nr, nc, bnnz, nnz, *bptr, *ptr, *row, *col, *index, *bindex;
  LIS_SCALAR* values;
  bool assembled;
} lis_matrix_t;

static void lis_solver_set_tolerances(void* context,
                                      real_t rel_tol,
                                      real_t abs_tol,
                                      real_t div_tol)

{
  lis_solver_t* solver = context;
  char command[129];
  snprintf(command, 128, "-conv_cond 0 -tol %g", rel_tol);
  lis_solver_set_option(command, solver->solver);
}

static void lis_solver_set_max_iterations(void* context,
                                          int max_iters)
{
  lis_solver_t* solver = context;
  char command[129];
  snprintf(command, 128, "-maxiter %d", max_iters);
  lis_solver_set_option(command, solver->solver);
}

static void lis_solver_set_operator(void* context,
                                    void* op)
{
  lis_solver_t* solver = context;
  lis_matrix_t* A = op;
  solver->op = A->A;
}

static void lis_solver_set_pc(void* context,
                              void* pc)
{
  lis_solver_t* solver = context;
  char* text = pc;
  lis_solver_set_option(text, solver->solver);
  polymec_free(text);
}

static bool lis_solver_solve(void* context,
                             void* b,
                             void* x,
                             real_t* res_norm,
                             int* num_iters)
{
  lis_solver_t* solver = context;
  LIS_VECTOR B = b;
  LIS_VECTOR X = x;
  polymec_suspend_fpe();
  LIS_INT err = lis_solve(solver->op, B, X, solver->solver);
  polymec_restore_fpe();
  bool solved = false;
  if (err == LIS_SUCCESS)
  {
    LIS_INT status;
    lis_solver_get_status(solver->solver, &status);
    solved = (status == LIS_SUCCESS);
    if (!solved)
    {
      if (status == LIS_BREAKDOWN)
        log_debug("krylov_solver_solve (LIS): division by zero");
      else if (status == LIS_OUT_OF_MEMORY)
        log_debug("krylov_solver_solve (LIS): out of memory");
      else if (status == LIS_MAXITER)
        log_debug("krylov_solver_solve (LIS): max # of iterations exceeded");
    }
  }
  lis_solver_get_residualnorm(solver->solver, res_norm);
  LIS_INT niters;
  lis_solver_get_iter(solver->solver, &niters);
  *num_iters = (int)niters;
  return solved;
}

static void lis_solver_dtor(void* context)
{
  lis_solver_t* solver = context;
  polymec_free(solver);
}

static krylov_solver_t* lis_factory_pcg_solver(void* context,
                                               MPI_Comm comm)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  lis_solver_set_option("-i cg", solver->solver);
  lis_solver_set_option("-p jacobi", solver->solver);
  lis_solver_set_option("-scale jacobi", solver->solver);
//  lis_solver_set_option("-print 2", solver->solver);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = lis_solver_set_tolerances,
                                 .set_max_iterations = lis_solver_set_max_iterations,
                                 .set_operator = lis_solver_set_operator,
                                 .set_preconditioner = lis_solver_set_pc,
                                 .solve = lis_solver_solve,
                                 .dtor = lis_solver_dtor};

  return krylov_solver_new("LIS PCG", solver, vtable);
}

static krylov_solver_t* lis_factory_gmres_solver(void* context,
                                                 MPI_Comm comm,
                                                 int krylov_dimension)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  char gmres[129];
  snprintf(gmres, 128, "-i gmres %d", krylov_dimension);
  lis_solver_set_option(gmres, solver->solver);
  lis_solver_set_option("-p jacobi", solver->solver);
  lis_solver_set_option("-scale jacobi", solver->solver);
//  lis_solver_set_option("-print 2", solver->solver);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = lis_solver_set_tolerances,
                                 .set_max_iterations = lis_solver_set_max_iterations,
                                 .set_operator = lis_solver_set_operator,
                                 .set_preconditioner = lis_solver_set_pc,
                                 .solve = lis_solver_solve,
                                 .dtor = lis_solver_dtor};
  return krylov_solver_new("LIS GMRES", solver, vtable);
}

static krylov_solver_t* lis_factory_bicgstab_solver(void* context,
                                                    MPI_Comm comm)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  lis_solver_set_option("-i bicgstab", solver->solver);
  lis_solver_set_option("-p jacobi", solver->solver);
  lis_solver_set_option("-scale jacobi", solver->solver);
//  lis_solver_set_option("-print 2", solver->solver);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = lis_solver_set_tolerances,
                                 .set_max_iterations = lis_solver_set_max_iterations,
                                 .set_operator = lis_solver_set_operator,
                                 .set_preconditioner = lis_solver_set_pc,
                                 .solve = lis_solver_solve,
                                 .dtor = lis_solver_dtor};
  return krylov_solver_new("LIS BiCGSTAB", solver, vtable);
}

static krylov_solver_t* lis_factory_special_solver(void* context,
                                                   MPI_Comm comm,
                                                   const char* solver_name,
                                                   string_string_unordered_map_t* options)
{
  POLYMEC_NOT_IMPLEMENTED
  return NULL;
}

static krylov_pc_t* lis_factory_pc(void* context,
                                   MPI_Comm comm,
                                   const char* pc_name,
                                   string_string_unordered_map_t* options)
{
  // We store the preconditioner as a text string.
  char* text = polymec_malloc(sizeof(char) * 129);
  snprintf(text, 128, "-p %s", pc_name);

  // Set up the virtual table.
  krylov_pc_vtable vtable = {.dtor = NULL}; 
  return krylov_pc_new(pc_name, text, vtable);
}

static lis_matrix_t* matrix_new(MPI_Comm comm, int block_size)
{
  ASSERT(block_size > 0);
  lis_matrix_t* mat = polymec_malloc(sizeof(lis_matrix_t));
  mat->comm = comm;
  mat->block_size = (LIS_INT)block_size;
  mat->nr = mat->nc = mat->bnnz = mat->nnz = 0;
  mat->bptr = NULL;
  mat->ptr = NULL;
  mat->row = NULL; 
  mat->col = NULL; 
  mat->index = NULL;
  mat->bindex = NULL;
  mat->values = NULL;
  mat->assembled = false;
  return mat;
}

static void matrix_set(lis_matrix_t* mat)
{
  LIS_INT type;
  lis_matrix_get_type(mat->A, &type);
  LIS_INT status;
  if (type == LIS_MATRIX_CSR)
    status = lis_matrix_set_csr(mat->nnz, mat->ptr, mat->index, mat->values, mat->A);
  else if (type == LIS_MATRIX_BSR)
    status = lis_matrix_set_bsr(mat->block_size, mat->block_size, mat->bnnz, mat->bptr, mat->bindex, mat->values, mat->A);
  else
  {
    ASSERT(type == LIS_MATRIX_VBR);
    status = lis_matrix_set_vbr(mat->nnz, mat->nr, mat->nc, mat->bnnz, mat->row, mat->col, mat->ptr, mat->bptr, mat->bindex, mat->values, mat->A);
  }
  ASSERT(status == LIS_SUCCESS);
  mat->assembled = false;
}

static void ensure_matrix_assembly(lis_matrix_t* mat)
{
  if (!mat->assembled)
  {
    lis_matrix_assemble(mat->A);
    mat->assembled = true;
  }
}

static void* matrix_clone(void* context)
{
  lis_matrix_t* mat = context;
  lis_matrix_t* clone = matrix_new(mat->comm, mat->block_size);
  lis_matrix_create(mat->comm, &clone->A);
  LIS_INT N_local, N_global;
  lis_matrix_get_size(mat->A, &N_local, &N_global);
  lis_matrix_set_size(clone->A, N_local, 0);
  LIS_INT type;
  lis_matrix_get_type(mat->A, &type);
  lis_matrix_set_type(clone->A, type);
  lis_matrix_set_destroyflag(clone->A, LIS_FALSE);
  clone->bnnz = mat->bnnz;
  clone->nnz = mat->nnz;
  clone->nr = mat->nr;
  clone->nc = mat->nc;
  clone->assembled = false;
  if (mat->bptr != NULL)
  {
    clone->bptr = polymec_malloc(sizeof(LIS_INT) * (clone->nr+1));
    memcpy(clone->bptr, mat->bptr, sizeof(LIS_INT) * (clone->nr+1));
  }
  if (mat->ptr != NULL)
  {
    int ptr_size = (type == LIS_MATRIX_CSR) ? (N_local+1) : (clone->bnnz+1);
    clone->ptr = polymec_malloc(sizeof(LIS_INT) * ptr_size);
    memcpy(clone->ptr, mat->ptr, sizeof(LIS_INT) * ptr_size);
  }
  if (mat->row != NULL)
  {
    clone->row = polymec_malloc(sizeof(LIS_INT) * (clone->nr+1));
    memcpy(clone->row, mat->row, sizeof(LIS_INT) * (clone->nr+1));
  }
  if (mat->col != NULL)
  {
    clone->col = polymec_malloc(sizeof(LIS_INT) * (clone->nc+1));
    memcpy(clone->col, mat->col, sizeof(LIS_INT) * (clone->nc+1));
  }
  if (mat->bindex != NULL)
  {
    clone->bindex = polymec_malloc(sizeof(LIS_INT) * clone->bnnz);
    memcpy(clone->bindex, mat->bindex, sizeof(LIS_INT) * clone->bnnz);
  }
  if (mat->index != NULL)
  {
    clone->index = polymec_malloc(sizeof(LIS_INT) * clone->nnz);
    memcpy(clone->index, mat->index, sizeof(LIS_INT) * clone->nnz);
  }
  if (mat->values != NULL)
  {
    clone->values = polymec_malloc(sizeof(LIS_SCALAR) * clone->nnz);
    memcpy(clone->values, mat->values, sizeof(LIS_SCALAR) * clone->nnz);
  }
  matrix_set(clone);
  ensure_matrix_assembly(clone);
  return clone;
}

static void matrix_zero(void* context)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);
  memset(mat->values, 0, sizeof(LIS_SCALAR) * mat->nnz);
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_scale(void* context, real_t scale_factor)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);
  for (int i = 0; i < mat->nnz; ++i)
    mat->values[i] *= scale_factor;
  matrix_set(mat);

  ensure_matrix_assembly(mat);
}

static void matrix_add_identity(void* context, real_t scale_factor)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);
  lis_matrix_shift_diagonal(mat->A, (LIS_SCALAR)scale_factor);
}

static void matrix_set_diagonal_csr(void* context, void* D)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  for (LIS_INT row = start; row < end; ++row)
  {
    LIS_SCALAR Di;
    lis_vector_get_value(d, row, &Di);
    LIS_INT row_offset = mat->ptr[row-start];
    LIS_INT ncol = mat->ptr[row+1-start] - row_offset;
    LIS_INT* col_start = &mat->index[row_offset];
    LIS_INT* col_ptr = (LIS_INT*)index_bsearch((index_t*)col_start, (int)ncol, (index_t)row);
    LIS_INT k = row_offset + (col_ptr - col_start);
    mat->values[k] = Di;
  }
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_add_diagonal_csr(void* context, void* D)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  for (LIS_INT row = start; row < end; ++row)
  {
    LIS_SCALAR Di;
    lis_vector_get_value(d, row, &Di);
    LIS_INT row_offset = mat->ptr[row-start];
    LIS_INT ncol = mat->ptr[row+1-start] - row_offset;
    LIS_INT* col_start = &mat->index[row_offset];
    LIS_INT* col_ptr = (LIS_INT*)index_bsearch((index_t*)col_start, (int)ncol, (index_t)row);
    LIS_INT k = row_offset + (col_ptr - col_start);
    mat->values[k] += Di;
  }
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_set_diagonal_bsr(void* context, void* D)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  LIS_INT bs = mat->block_size;
  start /= bs, end /= bs;
  for (LIS_INT i = start; i < end; ++i)
  {
    LIS_INT bk = mat->bptr[i];
    for (int j = 0; j < bs; ++j)
    {
      LIS_SCALAR Dij;
      lis_vector_get_value(d, bs*i+j, &Dij);
      mat->values[bk*bs*bs+bs*j+j] = Dij;
    }
  }
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_add_diagonal_bsr(void* context, void* D)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  int bs = mat->block_size;
  start /= bs, end /= bs;
  for (LIS_INT i = start; i < end; ++i)
  {
    LIS_INT bk = mat->bptr[i];
    for (int j = 0; j < bs; ++j)
    {
      LIS_SCALAR Dij;
      lis_vector_get_value(d, bs*i+j, &Dij);
      mat->values[bk*bs*bs+bs*j+j] += Dij;
    }
  }
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_set_diagonal_vbr(void* context, void* D)
{
  POLYMEC_NOT_IMPLEMENTED;
#if 0
  lis_matrix_t* mat = context;
  lis_matrix_unset(mat->A);
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  for (LIS_INT i = 0; i < mat->nr; ++i)
  {
    LIS_INT r = mat->row[i], r1 = mat->row[i+1];
    LIS_INT c = mat->col[r], c1 = mat->col[i+1];
    LIS_SCALAR Di;
    lis_vector_get_value(d, i, &Di);
    LIS_INT k = mat->ptr[i];
    mat->values[k] = Di;
  }
  matrix_set(mat);
  lis_matrix_assemble(mat->A);
#endif
}

static void matrix_add_diagonal_vbr(void* context, void* D)
{
  POLYMEC_NOT_IMPLEMENTED;
#if 0
  lis_matrix_t* mat = context;
  lis_matrix_unset(mat->A);
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  for (LIS_INT i = 0; i < mat->nr; ++i)
  {
    LIS_INT r = mat->row[i], r1 = mat->row[i+1];
    LIS_INT c = mat->col[r], c1 = mat->col[i+1];
    LIS_SCALAR Di;
    lis_vector_get_value(d, i, &Di);
    LIS_INT k = mat->ptr[i];
    mat->values[k] += Di;
  }
  matrix_set(mat);
  lis_matrix_assemble(mat->A);
#endif
}

static void matrix_set_values_csr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  lis_matrix_t* mat = context;
  lis_matrix_unset(mat->A);

  LIS_INT N_local, N_global;
  lis_matrix_get_size(mat->A, &N_local, &N_global);
  LIS_INT start, finish;
  lis_matrix_get_range(mat->A, &start, &finish);

  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    index_t row = rows[i];
    index_t num_cols = num_columns[i];

    if ((row < (index_t)start) || (row >= (index_t)finish)) continue; // skip off-proc rows.
    LIS_INT row_offset = mat->ptr[row-start];
    LIS_INT* col_start = &mat->index[row_offset];
    LIS_INT ncol = mat->ptr[row+1-start] - row_offset;

    for (index_t j = 0; j < num_cols; ++j, ++k)
    {
      // The columns are stored in our matrix in sorted order, 
      // so we can use a binary search to find the location of this one.
      index_t column = columns[k];
      LIS_INT* col_ptr = (LIS_INT*)index_bsearch((index_t*)col_start, (int)ncol, column);
      if (col_ptr == NULL) continue; // Column not found -- skip it.
      mat->values[row_offset+(col_ptr-col_start)] = values[k];
    }
  }
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_add_values_csr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  lis_matrix_t* mat = context;
  lis_matrix_unset(mat->A);

  LIS_INT N_local, N_global;
  lis_matrix_get_size(mat->A, &N_local, &N_global);
  LIS_INT start, finish;
  lis_matrix_get_range(mat->A, &start, &finish);

  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    index_t row = rows[i];
    index_t num_cols = num_columns[i];

    if ((row < (index_t)start) || (row >= (index_t)finish)) continue; // skip off-proc rows.
    LIS_INT row_offset = mat->ptr[row-start];
    LIS_INT* col_start = &mat->index[row_offset];
    LIS_INT ncol = mat->ptr[row+1-start] - row_offset;

    for (index_t j = 0; j < num_cols; ++j, ++k)
    {
      // The columns are stored in our matrix in sorted order, 
      // so we can use a binary search to find the location of this one.
      index_t column = columns[k];
      LIS_INT* col_ptr = (LIS_INT*)index_bsearch((index_t*)col_start, (int)ncol, column);
      if (col_ptr == NULL) continue; // Column not found -- skip it.
      mat->values[row_offset+(col_ptr-col_start)] += values[k];
    }
  }
  matrix_set(mat);
  ensure_matrix_assembly(mat);
}

static void matrix_get_values_csr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  LIS_INT N_local, N_global;
  lis_matrix_get_size(mat->A, &N_local, &N_global);
  LIS_INT start, finish;
  lis_matrix_get_range(mat->A, &start, &finish);

  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    index_t row = rows[i];
    index_t num_cols = num_columns[i];

    if ((row < (index_t)start) || (row >= (index_t)finish))
    {
      // All columns in this row are set to zero.
      memset(&values[k], 0, sizeof(real_t) * num_cols);
      continue;
    }

    LIS_INT row_offset = mat->ptr[row-start];
    LIS_INT* col_start = &mat->index[row_offset];
    LIS_INT ncol = mat->ptr[row+1-start] - row_offset;

    for (index_t j = 0; j < num_cols; ++j, ++k)
    {
      // The columns are stored in our matrix in sorted order, 
      // so we can use a binary search to find the location of this one.
      index_t column = columns[k];
      LIS_INT* col_ptr = (LIS_INT*)index_bsearch((index_t*)col_start, (int)ncol, column);
      if (col_ptr == NULL) // Column not found 
        values[k] = 0.0;
      else
        values[k] = mat->values[row_offset+(col_ptr-col_start)];
    }
  }
}

static void matrix_set_values_bsr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);

  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    int row = rows[i];
    for (int j = 0; j < num_columns[i]; ++j, ++k)
    {
      LIS_INT status = lis_matrix_set_value(LIS_INS_VALUE, row, columns[k], values[k], mat->A);
      ASSERT(status == LIS_SUCCESS);
    }
  }
  matrix_set(mat);
}

static void matrix_add_values_bsr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  lis_matrix_t* mat = context;
  ensure_matrix_assembly(mat);

  lis_matrix_unset(mat->A);

  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    int row = rows[i];
    for (int j = 0; j < num_columns[i]; ++j, ++k)
      lis_matrix_set_value(LIS_ADD_VALUE, row, columns[k], values[k], mat->A);
  }
  matrix_set(mat);
}

static void matrix_get_values_bsr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  POLYMEC_NOT_IMPLEMENTED
}

static void matrix_set_values_vbr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  POLYMEC_NOT_IMPLEMENTED
}

static void matrix_add_values_vbr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  POLYMEC_NOT_IMPLEMENTED
}

static void matrix_get_values_vbr(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  POLYMEC_NOT_IMPLEMENTED
}

static void matrix_dtor(void* context)
{
  lis_matrix_t* mat = context;
  lis_matrix_unset(mat->A);
  if (mat->bptr != NULL)
    polymec_free(mat->bptr);
  if (mat->ptr != NULL)
    polymec_free(mat->ptr);
  if (mat->row != NULL)
    polymec_free(mat->row);
  if (mat->col != NULL)
    polymec_free(mat->col);
  if (mat->bindex != NULL)
    polymec_free(mat->bindex);
  if (mat->index != NULL)
    polymec_free(mat->index);
  if (mat->values != NULL)
    polymec_free(mat->values);
  lis_matrix_destroy(mat->A);
  polymec_free(mat);
}

static krylov_matrix_t* lis_factory_matrix(void* context,
                                           matrix_sparsity_t* sparsity)
{
  MPI_Comm comm = matrix_sparsity_comm(sparsity);
  index_t N_local = matrix_sparsity_num_local_rows(sparsity);
  index_t N_global = matrix_sparsity_num_global_rows(sparsity);

  lis_matrix_t* mat = matrix_new(comm, 1);
  lis_matrix_create(comm, &mat->A);
  LIS_INT err = lis_matrix_set_size(mat->A, N_local, 0);
  if (err != LIS_SUCCESS)
    polymec_error("lis_factory_matrix: failed to create an %d x %d matrix.", N_global, N_global);
  lis_matrix_set_type(mat->A, LIS_MATRIX_CSR);
#ifndef NDEBUG
  {
    LIS_INT Nl, Ng;
    lis_matrix_get_size(mat->A, &Nl, &Ng);
    ASSERT(Nl == N_local);
    ASSERT(Ng == N_global);
  }
#endif

  // We manage all of these arrays ourself.
  lis_matrix_set_destroyflag(mat->A, LIS_FALSE);
  mat->nnz = matrix_sparsity_num_nonzeros(sparsity);
  mat->ptr = polymec_malloc(sizeof(LIS_INT) * (N_local+1));
  mat->index = polymec_malloc(sizeof(LIS_INT) * mat->nnz);
  mat->values = polymec_malloc(sizeof(LIS_SCALAR) * mat->nnz);

  index_t k = 0, row;
  int rpos = 0, r = 0;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    // Here's where the row begins.
    mat->ptr[r] = k;

    // We store the rows in strictly ascending column order, with
    // no exception for the diagonal.
    int cpos = 0;
    index_t col;
    while (matrix_sparsity_next_column(sparsity, row, &cpos, &col))
    {
      mat->index[k] = col;
      ++k;
    }

    // Sort this row.
    size_t nc = k - mat->ptr[r];
    index_qsort((index_t*)(&mat->index[mat->ptr[r]]), nc);
    ++r;
  }
  ASSERT(k == mat->nnz);
  mat->ptr[N_local] = mat->nnz;
  memset(mat->values, 0, sizeof(LIS_SCALAR) * mat->nnz);
  matrix_set(mat);
  lis_matrix_assemble(mat->A);
  mat->assembled = true;

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = matrix_clone,
                                 .zero = matrix_zero,
                                 .scale = matrix_scale,
                                 .add_identity = matrix_add_identity,
                                 .add_diagonal = matrix_add_diagonal_csr,
                                 .set_diagonal = matrix_set_diagonal_csr,
                                 .set_values = matrix_set_values_csr,
                                 .add_values = matrix_add_values_csr,
                                 .get_values = matrix_get_values_csr,
                                 .dtor = matrix_dtor};
  return krylov_matrix_new(mat, vtable, comm, N_local, N_global);
}

static krylov_matrix_t* lis_factory_block_matrix(void* context,
                                                 matrix_sparsity_t* sparsity,
                                                 int block_size)
{
  MPI_Comm comm = matrix_sparsity_comm(sparsity);
  index_t N_local = matrix_sparsity_num_local_rows(sparsity);
  index_t N_global = matrix_sparsity_num_global_rows(sparsity);

  // Initialize a Block Sparse Row (BSR) matrix with the given block size.
  lis_matrix_t* mat = matrix_new(comm, block_size);
  lis_matrix_create(comm, &mat->A);
  lis_matrix_set_size(mat->A, N_local, 0);
#if 0

  // We manage all of the arrays ourself.
  lis_matrix_set_destroyflag(mat->A, LIS_FALSE);
  int* adj = adj_graph_adjacency(sparsity);
  mat->bnnz = adj[N_local];
  mat->nnz = block_size * block_size * mat->bnnz;
  mat->nr = N_local, mat->nc = N_global;
  mat->bptr = polymec_malloc(sizeof(LIS_INT) * (mat->nr+1));
  mat->bindex = polymec_malloc(sizeof(LIS_INT) * mat->bnnz);
  mat->values = polymec_malloc(sizeof(LIS_SCALAR) * mat->nnz);
  int boffset = 0;
  for (int i = 0; i < mat->nr; ++i)
  {
    mat->bptr[i] = i;
    mat->bindex[boffset] = i;
    ++boffset;

    int ne = adj_graph_num_edges(sparsity, i);
    int* edges = adj_graph_edges(sparsity, i);
    for (int j = 0; j < ne; ++j, ++boffset)
    {
      int jj = edges[j];
      mat->bindex[boffset] = jj;
      memset(&mat->values[block_size*block_size*boffset], 0, sizeof(LIS_SCALAR) * block_size * block_size);
    }
  }
  mat->bptr[mat->nr] = mat->nr;
  matrix_set(mat);
  lis_matrix_assemble(mat->A);
  mat->assembled = true;

#endif
  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = matrix_clone,
                                 .zero = matrix_zero,
                                 .scale = matrix_scale,
                                 .add_identity = matrix_add_identity,
                                 .add_diagonal = matrix_add_diagonal_bsr,
                                 .set_diagonal = matrix_set_diagonal_bsr,
                                 .set_values = matrix_set_values_bsr,
                                 .add_values = matrix_add_values_bsr,
                                 .get_values = matrix_get_values_bsr,
                                 .dtor = matrix_dtor};
  return krylov_matrix_new(mat, vtable, comm, N_local, N_global);
}

static krylov_matrix_t* lis_factory_var_block_matrix(void* context,
                                                     matrix_sparsity_t* sparsity,
                                                     int* block_sizes)
{
  MPI_Comm comm = matrix_sparsity_comm(sparsity);
  index_t N_local = matrix_sparsity_num_local_rows(sparsity);
  index_t N_global = matrix_sparsity_num_global_rows(sparsity);

  lis_matrix_t* mat = matrix_new(comm, 1);
  lis_matrix_create(comm, &mat->A);
  lis_matrix_set_size(mat->A, N_local, 0);

#if 0
  // We manage all of the arrays ourself.
  lis_matrix_set_destroyflag(mat->A, LIS_FALSE);
  mat->nr = N_local; 
  mat->nc = N_global;
  int* adj = adj_graph_adjacency(sparsity);
  mat->bnnz = adj[N_local];
  mat->row = polymec_malloc(sizeof(LIS_INT) * (mat->nr+1));
  mat->col = polymec_malloc(sizeof(LIS_INT) * (mat->nc+1));
  mat->ptr = polymec_malloc(sizeof(LIS_INT) * (mat->bnnz+1));
  mat->bindex = polymec_malloc(sizeof(LIS_INT) * mat->bnnz);
  mat->bptr = polymec_malloc(sizeof(LIS_INT) * (mat->nr+1));
  LIS_INT boffset = 0, row_offset = 0;
  mat->nnz = 0;
  for (int i = 0; i < N_local; ++i)
  {
    int block_size = block_sizes[i];
    mat->row[i] = row_offset;
    mat->col[i] = row_offset;
    row_offset += block_size;
    mat->bptr[i] = i;
    mat->nnz += block_size * block_size;

    int ne = adj_graph_num_edges(sparsity, i);
    int* edges = adj_graph_edges(sparsity, i);
    for (int j = 0; j < ne; ++j, ++boffset)
    {
      int jj = edges[j];
      mat->bindex[boffset] = jj;
      mat->ptr[boffset] = mat->nnz;
      mat->nnz += block_size * block_size;
    }
  }
  mat->ptr[mat->bnnz] = mat->nnz;
  mat->values = polymec_malloc(sizeof(LIS_SCALAR) * mat->nnz);
  memset(mat->values, 0, sizeof(LIS_SCALAR) * mat->nnz);
  matrix_set(mat);
  lis_matrix_assemble(mat->A);
  mat->assembled = true;
#endif

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = matrix_clone,
                                 .zero = matrix_zero,
                                 .scale = matrix_scale,
                                 .add_identity = matrix_add_identity,
                                 .add_diagonal = matrix_add_diagonal_vbr,
                                 .set_diagonal = matrix_set_diagonal_vbr,
                                 .set_values = matrix_set_values_vbr,
                                 .add_values = matrix_add_values_vbr,
                                 .get_values = matrix_get_values_vbr,
                                 .dtor = matrix_dtor};
  return krylov_matrix_new(mat, vtable, comm, N_local, N_global);
}

static void* vector_clone(void* context)
{
  LIS_VECTOR v = context;
  LIS_VECTOR clone;
  lis_vector_duplicate(v, &clone);
  lis_vector_copy(v, clone);
  return clone;
}

static void vector_zero(void* context)
{
  LIS_VECTOR v = context;
  lis_vector_set_all(0.0, v);
}

static void vector_set_value(void* context, real_t value)
{
  LIS_VECTOR v = context;
  lis_vector_set_all((LIS_SCALAR)value, v);
}

static void vec_scale(void* context, real_t scale_factor)
{
  LIS_VECTOR v = context;
  lis_vector_scale((LIS_SCALAR)scale_factor, v);
}

static void vector_set_values(void* context, index_t num_values,
                              index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  LIS_INT ind[num_values];
  for (index_t i = 0; i < num_values; ++i)
    ind[i] = (LIS_INT)indices[i];
  LIS_INT status = lis_vector_set_values(LIS_INS_VALUE, (LIS_INT)num_values, 
                                         ind, values, v);
  ASSERT(status == LIS_SUCCESS);
}

static void vector_add_values(void* context, index_t num_values,
                              index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  LIS_INT ind[num_values];
  for (index_t i = 0; i < num_values; ++i)
    ind[i] = (LIS_INT)indices[i];
  LIS_INT status = lis_vector_set_values(LIS_ADD_VALUE, (LIS_INT)num_values, 
                                         ind, values, v);
  ASSERT(status == LIS_SUCCESS);
}

static void vector_get_values(void* context, index_t num_values,
                              index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  for (index_t i = 0; i < num_values; ++i)
    lis_vector_get_values(v, (LIS_INT)indices[i], 1, &values[i]);
}

static real_t vector_norm(void* context, int p)
{
  LIS_VECTOR v = context;
  real_t norm;
  if (p == 0)
    lis_vector_nrmi(v, &norm);
  else if (p == 1)
    lis_vector_nrm1(v, &norm);
  else // (p == 2)
    lis_vector_nrm2(v, &norm);
  return norm;
}

static void vector_dtor(void* context)
{
  LIS_VECTOR v = context;
  lis_vector_destroy(v);
}

static krylov_vector_t* lis_factory_vector(void* context,
                                           MPI_Comm comm,
                                           index_t* row_dist)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  index_t N_local = row_dist[rank+1] - row_dist[rank];
  index_t N_global = row_dist[nprocs];

  LIS_VECTOR v;
  lis_vector_create(comm, &v);
  lis_vector_set_size(v, N_local, 0);
  lis_vector_set_all(0.0, v);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = vector_clone,
                                 .zero = vector_zero,
                                 .set_value = vector_set_value,
                                 .scale = vec_scale,
                                 .set_values = vector_set_values,
                                 .add_values = vector_add_values,
                                 .get_values = vector_get_values,
                                 .norm = vector_norm,
                                 .dtor = vector_dtor};
  return krylov_vector_new(v, vtable, N_local, N_global);
}

static void lis_factory_dtor(void* context)
{
  lis_factory_t* factory = context;
  polymec_free(factory);
}

krylov_factory_t* lis_krylov_factory()
{
  _Static_assert(sizeof(LIS_INT) == sizeof(index_t), "sizeof(LIS_INT) != sizeof(index_t)");
  _Static_assert(sizeof(LIS_REAL) == sizeof(real_t), "LIS_REAL != real_t");

  lis_factory_t* factory = polymec_malloc(sizeof(lis_factory_t));

  // Construct the factory.
  krylov_factory_vtable vtable = {.pcg_solver = lis_factory_pcg_solver,
                                  .gmres_solver = lis_factory_gmres_solver,
                                  .bicgstab_solver = lis_factory_bicgstab_solver,
                                  .special_solver = lis_factory_special_solver,
                                  .preconditioner = lis_factory_pc,
                                  .matrix = lis_factory_matrix,
                                  .block_matrix = lis_factory_block_matrix,
                                  .var_block_matrix = lis_factory_var_block_matrix,
                                  .vector = lis_factory_vector,
                                  .dtor = lis_factory_dtor};
  return krylov_factory_new("LIS", factory, vtable);
}

