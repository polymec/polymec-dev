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

#include "core/linear_algebra.h"
#include "core/sundials_helpers.h"
#include "integrators/block_jacobi_preconditioner.h"

typedef struct 
{
  int (*F)(void* context, real_t t, real_t* x, real_t* F);
  void (*communicate)(void* context, real_t t, real_t* x);
  void* context;

  int num_block_rows;
  int block_size; 

  // Work vectors.
  int num_work_vectors;
  real_t** work;

} block_jacobi_preconditioner_t;

// Block-diagonal matrix.
typedef struct
{
  int num_block_rows;
  int block_size;
  real_t* coeffs;
} bd_mat_t;

static void bd_scale_and_shift(void* context, real_t gamma)
{
  bd_mat_t* mat = context;
  for (int i = 0; i < mat->num_block_rows*mat->block_size; ++i)
    mat->coeffs[i] = 1.0 + gamma * mat->coeffs[i];
}

static real_t bd_coeff(void* context, int i, int j)
{
  bd_mat_t* mat = context;
  int bs = mat->block_size;
  if (abs(j - i) >= bs)
    return 0.0;
  int block_row = i/bs;
  int r = i % bs;
  int c = j - block_row*bs;
  real_t* A = &mat->coeffs[block_row * bs * bs];
  return A[bs * r + c];
}

static void bd_dtor(void* context)
{
  bd_mat_t* mat = context;
  free(mat->coeffs);
  free(mat);
}

static preconditioner_matrix_t* block_jacobi_preconditioner_matrix(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  preconditioner_matrix_vtable vtable = {.scale_and_shift = bd_scale_and_shift,
                                         .coeff = bd_coeff,
                                         .dtor = bd_dtor};
  int bs = precond->block_size;
  int num_rows = precond->num_block_rows * bs;
  bd_mat_t* mat = malloc(sizeof(bd_mat_t));
  mat->num_block_rows = precond->num_block_rows;
  mat->block_size = bs;
  mat->coeffs = malloc(sizeof(real_t) * mat->num_block_rows * bs * bs);
  memset(mat->coeffs, 0, sizeof(real_t) * mat->num_block_rows * bs * bs);
  return preconditioner_matrix_new("Block-diagonal", mat, vtable, num_rows);
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

static void insert_Jv_into_bd_mat(int num_rows, 
                                  int color, 
                                  real_t* Jv, 
                                  bd_mat_t* J)
{
  for (int i = color; i < num_rows; i += J->block_size)
  {
    int block_row = i / J->block_size;
    for (int j = block_row * J->block_size; j < (block_row+1)*J->block_size; ++j)
    {
      int c = j % J->block_size;
      J->coeffs[block_row*J->block_size+c] = Jv[j];
    }
  }
}

static void block_jacobi_preconditioner_compute_jacobian(void* context, real_t t, real_t* x, preconditioner_matrix_t* mat)
{
  block_jacobi_preconditioner_t* precond = context;
  bd_mat_t* A = preconditioner_matrix_context(mat);
  real_t** work = precond->work;

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int num_rows = A->num_block_rows * A->block_size;
  real_t* Jv = malloc(sizeof(real_t) * num_rows);
  int num_colors = A->block_size;
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * num_rows);
    for (int i = c; i < num_rows; i += A->block_size)
      work[0][i] = 1.0;

    // We evaluate F(x) and place it into work[1].
    if (precond->communicate != NULL)
      precond->communicate(precond->context, t, x);
    precond->F(precond->context, t, x, work[1]);

    // Now evaluate the matrix-vector product.
    memset(Jv, 0, sizeof(real_t) * num_rows);
    finite_diff_Jv(precond->F, precond->communicate, precond->context, x, t, num_rows, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    insert_Jv_into_bd_mat(num_rows, c, Jv, A);
  }
  free(Jv);
}

static bool block_jacobi_preconditioner_solve(void* context, preconditioner_matrix_t* A, real_t* B)
{
  block_jacobi_preconditioner_t* precond = context;
  int bs = precond->block_size;
  bd_mat_t* mat = preconditioner_matrix_context(A);

  for (int i = 0; i < precond->num_block_rows; ++i)
  {
    // Copy the block for this row into place.
    real_t Aij[bs*bs], bi[bs];
    memcpy(Aij, &mat->coeffs[i*bs*bs], sizeof(real_t)*bs*bs);
    memcpy(bi, &B[i*bs], sizeof(real_t)*bs);

    // Solve the linear system.
    int one = 1, ipiv[bs], info;
    dgesv(&bs, &one, Aij, &bs, ipiv, bi, &bs, &info);
    ASSERT(info == 0);

    // Copy the solution into place.
    memcpy(&B[i*bs], bi, sizeof(real_t)*bs);
  }

  return true;
}

static void block_jacobi_preconditioner_dtor(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  for (int i = 0; i < precond->num_work_vectors; ++i)
    free(precond->work[i]);
  free(precond->work);
  free(precond);
}

preconditioner_t* block_jacobi_preconditioner_new(void* context,
                                                  int (*residual_func)(void* context, real_t t, real_t* x, real_t* F),
                                                  void (*communication_func)(void* context, real_t t, real_t* x),
                                                  int num_block_rows,
                                                  int block_size)
{
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);

  block_jacobi_preconditioner_t* precond = malloc(sizeof(block_jacobi_preconditioner_t));
  precond->F = residual_func;
  precond->communicate = communication_func;
  precond->context = context;
  precond->num_block_rows = num_block_rows;
  precond->block_size = block_size;

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = malloc(sizeof(real_t) * precond->num_block_rows * precond->block_size);

  preconditioner_vtable vtable = {.matrix = block_jacobi_preconditioner_matrix,
                                  .compute_jacobian = block_jacobi_preconditioner_compute_jacobian,
                                  .solve = block_jacobi_preconditioner_solve,
                                  .dtor = block_jacobi_preconditioner_dtor};
  return preconditioner_new("Block Jacobi preconditioner", precond, vtable);
}
