// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"
#include "solvers/bd_matrix.h"

struct bd_matrix_t
{
  size_t num_block_rows;
  int *D_offsets, *B_offsets; // For variable block sizes.
  real_t* D;
  int block_size; // -1 if variable, set to constant size if applicable.
};

bd_matrix_t* bd_matrix_new(size_t num_block_rows,
                           size_t block_size)
{
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);
  size_t block_sizes[num_block_rows];
  for (size_t i = 0; i < num_block_rows; ++i)
    block_sizes[i] = block_size;
  return var_bd_matrix_new(num_block_rows, block_sizes);
}

bd_matrix_t* var_bd_matrix_new(size_t num_block_rows,
                               size_t* block_sizes)
{
  ASSERT(num_block_rows > 0);

  bd_matrix_t* A = polymec_malloc(sizeof(bd_matrix_t));
  A->num_block_rows = num_block_rows;
  A->D_offsets = polymec_malloc(sizeof(int) * (num_block_rows+1));
  A->B_offsets = polymec_malloc(sizeof(int) * (num_block_rows+1));
  A->D_offsets[0] = A->B_offsets[0] = 0;
  bool constant_block_size = true;
  int bs0 = -1;
  for (int i = 0; i < (int)num_block_rows; ++i)
  {
    size_t bs = block_sizes[i];
    if (bs0 == -1)
      bs0 = (int)bs;
    else if (bs != bs0)
      constant_block_size = false;
    ASSERT(bs >= 1);
    A->D_offsets[i+1] = (int)(A->D_offsets[i] + bs*bs);
    A->B_offsets[i+1] = (int)(A->B_offsets[i] + bs);
  }
  if (constant_block_size)
    A->block_size = bs0;
  else
    A->block_size = -1;
  int N = A->D_offsets[A->num_block_rows];
  A->D = polymec_calloc(N, sizeof(real_t));
  return A;
}

bd_matrix_t* bd_matrix_clone(bd_matrix_t* matrix)
{
  bd_matrix_t* clone = polymec_malloc(sizeof(bd_matrix_t));
  clone->num_block_rows = matrix->num_block_rows;
  clone->D_offsets = polymec_malloc(sizeof(int) * (clone->num_block_rows+1));
  memcpy(clone->D_offsets, matrix->D_offsets, sizeof(int) * (clone->num_block_rows+1));
  clone->B_offsets = polymec_malloc(sizeof(int) * (clone->num_block_rows+1));
  memcpy(clone->B_offsets, matrix->B_offsets, sizeof(int) * (clone->num_block_rows+1));
  int N = clone->D_offsets[clone->num_block_rows];
  clone->D = polymec_malloc(sizeof(real_t) * N);
  memcpy(clone->D, matrix->D, sizeof(real_t) * N);
  clone->block_size = matrix->block_size;
  return clone;
}

void bd_matrix_free(bd_matrix_t* matrix)
{
  polymec_free(matrix->D);
  polymec_free(matrix->D_offsets);
  polymec_free(matrix->B_offsets);
  polymec_free(matrix);
}

size_t bd_matrix_num_block_rows(bd_matrix_t* matrix)
{
  return matrix->num_block_rows;
}

size_t bd_matrix_num_rows(bd_matrix_t* matrix)
{
  return matrix->B_offsets[matrix->num_block_rows];
}

void bd_matrix_zero(bd_matrix_t* matrix)
{
  memset(matrix->D, 0, sizeof(real_t) * matrix->D_offsets[matrix->num_block_rows]);
}

void bd_matrix_add_identity(bd_matrix_t* matrix, real_t scale_factor)
{
  for (int i = 0; i < matrix->num_block_rows; ++i)
  {
    int bs = matrix->B_offsets[i+1] - matrix->B_offsets[i];
    int offset = matrix->D_offsets[i];
    for (int j = 0; j < bs; ++j)
      matrix->D[offset+bs*j+j] += scale_factor;
  }
}

size_t bd_matrix_block_size(bd_matrix_t* matrix, int block_row)
{
  return (matrix->block_size != -1) ? (size_t)matrix->block_size 
                                    : (size_t)(matrix->B_offsets[block_row+1] - matrix->B_offsets[block_row]);
}

void bd_matrix_insert_block(bd_matrix_t* matrix, int block_row, real_t* block)
{
  size_t bs = bd_matrix_block_size(matrix, block_row);
  real_t* B = bd_matrix_block(matrix, block_row);
  memcpy(B, block, sizeof(real_t) * bs * bs);
}

void bd_matrix_add_block(bd_matrix_t* matrix, int block_row, real_t* block)
{
  size_t bs = bd_matrix_block_size(matrix, block_row);
  real_t* B = bd_matrix_block(matrix, block_row);
  for (int i = 0; i < bs*bs; ++i)
    B[i] += block[i];
}

void bd_matrix_matvec(bd_matrix_t* matrix, real_t* vector, real_t* product)
{
  for (int i = 0; i < (int)matrix->num_block_rows; ++i)
  {
    int bs = matrix->B_offsets[i+1] - matrix->B_offsets[i];
    real_t* Ai = &matrix->D[matrix->D_offsets[i]];
    real_t* xi = &vector[matrix->B_offsets[i]];
    real_t* Axi = &product[matrix->B_offsets[i]];
    char no_trans = 'N';
    real_t one = 1.0, zero = 0.0;
    int incx = 1;
    rgemv(&no_trans, &bs, &bs, &one, Ai, &bs, xi, &incx, &zero, Axi, &incx);
  }
}

void bd_matrix_fprintf(bd_matrix_t* matrix, FILE* stream)
{
  size_t N = matrix->num_block_rows;
  real_t* D = matrix->D;
  fprintf(stream, "\nBlock diagonal matrix (N = %d):\n", (int)N);
  for (int i = 0; i < (int)N; ++i)
  {
    fprintf(stream, "%d: [", i);
    int bs = matrix->B_offsets[i+1] - matrix->B_offsets[i];
    int offset = matrix->D_offsets[i];
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

real_t* bd_matrix_block(bd_matrix_t* matrix, int block_row)
{
  return &(matrix->D[matrix->D_offsets[block_row]]);
}

bool solve_bd_system(bd_matrix_t* A, real_t* B, real_t* x)
{
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
      if (reals_equal(Aij[bs*j+j], 0.0))
        Aij[bs*j+j] = epsilon;
    }

    // Solve the linear system.
    int one = 1, ipiv[bs], info;
    rgesv(&bs, &one, Aij, &bs, ipiv, bi, &bs, &info);
    success = (info == 0);

    if (success)
    {
      // Copy the solution into place.
      memcpy(&x[B_offset], bi, sizeof(real_t)*bs);
    }
    else
    {
      ASSERT(info > 0);
      log_debug("bd_matrix_solve: call to rgesv failed for block row %d.", i);
      log_debug("bd_matrix_solve: (U is singular.)");
      break;
    }
  }

  return success;
}

