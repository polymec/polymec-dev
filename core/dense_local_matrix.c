// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"
#include "core/dense_local_matrix.h"

typedef struct 
{
  int N; // Number of rows in matrix.

  // Matrix data.
  real_t* A;
} dlm_t;

static void* dlm_clone(void* context)
{
  dlm_t* mat = context;
  dlm_t* clone = polymec_malloc(sizeof(dlm_t));
  clone->N = mat->N;
  clone->A = polymec_malloc(sizeof(real_t) * mat->N * mat->N);
  memcpy(clone->A, mat->A, sizeof(real_t) * mat->N * mat->N);
  return clone;
}

static void dlm_zero(void* context)
{
  dlm_t* mat = context;
  memset(mat->A, 0, sizeof(real_t) * mat->N * mat->N);
}

static void dlm_add_identity(void* context, real_t scale_factor)
{
  dlm_t* mat = context;

  // Add in the diagonal values.
  for (int i = 0; i < mat->N; ++i)
    mat->A[mat->N*i + i] += scale_factor;
}

static int dlm_num_columns(void* context, int row)
{
  dlm_t* mat = context;
  return mat->N;
}

static void dlm_get_columns(void* context, int row, int* columns)
{
  dlm_t* mat = context;
  for (int j = 0; j < mat->N; ++j)
    columns[j] = mat->A[mat->N*j+row];
}

static void dlm_add_column_vector(void* context, 
                                  real_t scale_factor,
                                  int column,
                                  real_t* column_vector)
{
  dlm_t* mat = context;

  if (column >= mat->N) return;

  // Add in the row values.
  for (int i = 0; i < mat->N; ++i)
    mat->A[mat->N*column+i] += scale_factor * column_vector[i];
}

static void dlm_add_row_vector(void* context, 
                               real_t scale_factor,
                               int row,
                               real_t* row_vector)
{
  dlm_t* mat = context;

  if (row >= mat->N) return;

  // Add in the column values.
  for (int j = 0; j < mat->N; ++j)
    mat->A[mat->N*j+row] += scale_factor * row_vector[j];
}

static bool dlm_solve(void* context, real_t* B, real_t* x)
{
  dlm_t* mat = context;

  // Solve the linear system.
  int one = 1, ipiv[mat->N], info;
  memcpy(x, B, sizeof(real_t) * mat->N);
  rgesv(&mat->N, &one, mat->A, &mat->N, ipiv, x, &mat->N, &info);
  bool success = (info == 0);

  if (!success)
  {
    ASSERT(info > 0);
    log_debug("dlm_solve: call to dgesv failed.");
    log_debug("dlm_solve: (U is singular.)");
  }

  return success;
}

static void dlm_fprintf(void* context, FILE* stream)
{
  dlm_t* mat = context;
  int N = mat->N;
  fprintf(stream, "\nDense matrix (N = %d):\n", mat->N);
  matrix_fprintf(mat->A, N, N, stream);
}

static real_t dlm_value(void* context, int i, int j)
{
  dlm_t* mat = context;
  ASSERT(i < mat->N);
  ASSERT(j < mat->N);
  return mat->A[mat->N*j+i];
}

static void dlm_set_value(void* context, int i, int j, real_t value)
{
  dlm_t* mat = context;
  ASSERT(i < mat->N);
  ASSERT(j < mat->N);
  mat->A[mat->N*j+i] = value;
}

static void dlm_get_diag(void* context, real_t* diag)
{
  dlm_t* mat = context;
  for (int i = 0; i < mat->N; ++i)
    diag[i] = mat->A[mat->N*i+i];
}

static void dlm_matvec(void* context, real_t* x, real_t* Ax)
{
  dlm_t* mat = context;
  char no_trans = 'N';
  real_t one = 1.0, zero = 0.0;
  int incx = 1;
  rgemv(&no_trans, &mat->N, &mat->N, &one, mat->A, &mat->N, x, &incx, &zero, Ax, &incx);
}

static void dlm_add_matrix(void* context, real_t scale_factor, void* B)
{
  dlm_t* Amat = context;
  dlm_t* Bmat = B;
  for (int i = 0; i < Amat->N*Amat->N; ++i)
    Amat->A[i] += scale_factor * Bmat->A[i];
}

static real_t dlm_norm(void* context, char n)
{
  dlm_t* mat = context;
  real_t work[mat->N];
  return rlange(&n, &mat->N, &mat->N, mat->A, &mat->N, work);
}

static void dlm_dtor(void* context)
{
  dlm_t* mat = context;
  if (mat->A != NULL)
    polymec_free(mat->A);
  polymec_free(mat);
}

local_matrix_t* dense_local_matrix_new(int N)
{
  ASSERT(N > 0);
  dlm_t* mat = polymec_malloc(sizeof(dlm_t));
  mat->N = N;
  mat->A = polymec_malloc(sizeof(real_t) * N * N);
  memset(mat->A, 0, sizeof(real_t) * N * N);

  char name[1024];
  snprintf(name, 1024, "Dense local matrix (N = %d)", mat->N);
  local_matrix_vtable vtable = {.clone = dlm_clone,
                                .dtor = dlm_dtor,
                                .zero = dlm_zero,
                                .num_columns = dlm_num_columns,
                                .get_columns = dlm_get_columns,
                                .add_identity = dlm_add_identity,
                                .add_column_vector = dlm_add_column_vector,
                                .add_row_vector = dlm_add_row_vector,
                                .solve = dlm_solve,
                                .fprintf = dlm_fprintf,
                                .value = dlm_value,
                                .set_value = dlm_set_value,
                                .get_diag = dlm_get_diag,
                                .matvec = dlm_matvec,
                                .add_matrix = dlm_add_matrix,
                                .norm = dlm_norm};
  return local_matrix_new(name, mat, vtable, N);
}

