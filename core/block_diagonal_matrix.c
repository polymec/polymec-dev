// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/block_diagonal_matrix.h"
#include "core/linear_algebra.h"

typedef struct 
{
  int num_block_rows;
  int *D_offsets, *B_offsets; // For variable block sizes.
  real_t* D;
  int block_size; // -1 if variable, set to constant size if applicable.
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

  // We have to find the right block column.
  int block_col = 0;
  while ((A->B_offsets[block_col] < i) && (block_col < A->num_block_rows))
    ++block_col;

  if (block_col < A->num_block_rows)
  {
    int bs = A->B_offsets[block_col+1] - A->B_offsets[block_col];
    int c = i % bs;
    for (int r = 0; r < bs; ++r)
    {
      int j = A->B_offsets[block_col] + r;
      D[A->D_offsets[block_col] + c*bs + r] += scale_factor * column_vector[j];
    }
  }
}

static void bdm_add_column_vector_constant_bs(void* context,
                                              real_t scale_factor,
                                              int column,
                                              real_t* column_vector)
{
  bdm_t* A = context;
  ASSERT(A->block_size > 0);
  real_t* D = A->D;

  int i = column;
  if (i < A->B_offsets[A->num_block_rows])
  {
    int bs = A->block_size;
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
      log_debug("bdm_solve: (U is singular.)");
      break;
    }
  }

  return success;
}

static void bdm_fprintf(void* context, FILE* stream)
{
  bdm_t* A = context;
  int N = A->num_block_rows;
  real_t* D = A->D;
  fprintf(stream, "\nBlock diagonal matrix (N = %d):\n", N);
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

static real_t bdm_value(void* context, int i, int j)
{
  bdm_t* A = context;
  real_t* D = A->D;

  // We have to find the right block column.
  int block_col = 0;
  while ((A->B_offsets[block_col+1] < i) && (block_col < A->num_block_rows))
    ++block_col;

  if (block_col < A->num_block_rows)
  {
    int bs = A->B_offsets[block_col+1] - A->B_offsets[block_col];
    int r = i % bs;
    if ((j >= block_col*bs) && (j < (block_col+1)*bs))
    {
      int c = j % bs;
      return D[A->D_offsets[block_col] + c*bs + r];
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

static real_t bdm_value_constant_bs(void* context, int i, int j)
{
  bdm_t* A = context;
  ASSERT(A->block_size > 0);
  real_t* D = A->D;

  if (i < A->B_offsets[A->num_block_rows])
  {
    int bs = A->block_size;
    int block_col = i / bs;
    int r = i % bs;
    if ((j >= block_col*bs) && (j < (block_col+1)*bs))
    {
      int c = j % bs;
      return D[A->D_offsets[block_col] + c*bs + r];
    }
    else
      return 0.0;
  }
  else
    return 0.0;
}

static void bdm_set_value(void* context, int i, int j, real_t value)
{
  bdm_t* A = context;
  real_t* D = A->D;

  // We have to find the right block column.
  int block_col = 0;
  while ((A->B_offsets[block_col+1] < i) && (block_col < A->num_block_rows))
    ++block_col;

  if (block_col < A->num_block_rows)
  {
    int bs = A->B_offsets[block_col+1] - A->B_offsets[block_col];
    int r = i % bs;
    if ((j >= block_col*bs) && (j < (block_col+1)*bs))
    {
      int c = j % bs;
      D[A->D_offsets[block_col] + c*bs + r] = value;
    }
  }
}

static void bdm_set_value_constant_bs(void* context, int i, int j, real_t value)
{
  bdm_t* A = context;
  ASSERT(A->block_size > 0);
  real_t* D = A->D;

  if (i < A->B_offsets[A->num_block_rows])
  {
    int bs = A->block_size;
    int block_col = i / bs;
    int r = i % bs;
    if ((j >= block_col*bs) && (j < (block_col+1)*bs))
    {
      int c = j % bs;
      D[A->D_offsets[block_col] + c*bs + r] = value;
    }
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
  bool constant_block_size = true;
  int bs0 = -1;
  for (int i = 0; i < num_block_rows; ++i)
  {
    int bs = block_sizes[i];
    if (bs0 == -1)
      bs0 = bs;
    else if (bs != bs0)
      constant_block_size = false;
    ASSERT(bs >= 1);
    A->D_offsets[i+1] = A->D_offsets[i] + bs*bs;
    A->B_offsets[i+1] = A->B_offsets[i] + bs;
  }
  if (constant_block_size)
    A->block_size = bs0;
  else
    A->block_size = -1;
  int N = A->D_offsets[A->num_block_rows];
  A->D = polymec_malloc(sizeof(real_t) * N);

  char name[1024];
  if (constant_block_size)
    snprintf(name, 1024, "Block diagonal matrix (bs = %d)", bs0);
  else
    snprintf(name, 1024, "Variable block diagonal matrix");
  local_matrix_vtable vtable = {.dtor = bdm_dtor,
                                .zero = bdm_zero,
                                .add_identity = bdm_add_identity,
                                .add_column_vector = bdm_add_column_vector,
                                .solve = bdm_solve,
                                .fprintf = bdm_fprintf,
                                .value = bdm_value,
                                .set_value = bdm_set_value};
  if (constant_block_size)
  {
    vtable.add_column_vector = bdm_add_column_vector_constant_bs;
    vtable.value = bdm_value_constant_bs;
    vtable.set_value = bdm_set_value_constant_bs;
  }
  return local_matrix_new(name, A, vtable);
}

