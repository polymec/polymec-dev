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
#include "core/block_jacobi_preconditioner.h"

typedef struct 
{
  void* context;

  int num_block_rows;
  int block_size; 

  // (Block) diagonal.
  real_t* D;

  // Method of computing the diagonal.
  void (*compute_diagonal)(void* context, int block_size, real_t* D);

  // Destructor for context pointer.
  void (*dtor)(void* context);
} block_jacobi_preconditioner_t;

static void block_jacobi_preconditioner_setup(void* context)
{
  block_jacobi_preconditioner_t* precond = context;

  // Zero our preconditioner matrix.
  int n = precond->num_block_rows;
  int bs = precond->block_size;
  memset(precond->D, 0, sizeof(real_t) * n * bs * bs);

  // Compute the block diagonal.
  precond->compute_diagonal(precond->context, precond->block_size, precond->D);
}


static bool block_jacobi_preconditioner_solve(void* context, real_t* B)
{
  block_jacobi_preconditioner_t* precond = context;
  int bs = precond->block_size;
  real_t* D = precond->D;

  bool success = false;
  for (int i = 0; i < precond->num_block_rows; ++i)
  {
    // Copy the block for this row into place.
    real_t Aij[bs*bs], bi[bs];
    memcpy(Aij, &D[i*bs*bs], sizeof(real_t)*bs*bs);
    memcpy(bi, &B[i*bs], sizeof(real_t)*bs);

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
      memcpy(&B[i*bs], bi, sizeof(real_t)*bs);
    }
    else
    {
      ASSERT(info > 0);
      log_debug("block_jacobi_preconditioner_solve: call to dgesv failed for block row %d.", i);
      log_debug("(U is singular).", i);
      break;
    }
  }

  return success;
}

static void block_jacobi_preconditioner_fprintf(void* context, FILE* stream)
{
  block_jacobi_preconditioner_t* precond = context;
  int n = precond->num_block_rows;
  int bs = precond->block_size;
  fprintf(stream, "Block Jacobian preconditioner matrix: (%d rows):", n * bs);
  for (int i = 0; i < n; ++i)
  {
    fprintf(stream, "\nRows %6d - %6d: ", i*bs, (i+1)*bs - 1);
    matrix_fprintf(&precond->D[i * bs * bs], bs, bs, stream);
  }
  fprintf(stream, "\n");
}

static void block_jacobi_preconditioner_dtor(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  polymec_free(precond->D);
  if ((precond->dtor != NULL) && (precond->context != NULL))
    precond->dtor(precond->context);
  polymec_free(precond);
}

preconditioner_t* block_jacobi_preconditioner_new(void* context,
                                                  void (*compute_diagonal)(void* context, int block_size, real_t* D),
                                                  void (*dtor)(void* context),
                                                  int num_block_rows,
                                                  int block_size)
{
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);

  block_jacobi_preconditioner_t* precond = polymec_malloc(sizeof(block_jacobi_preconditioner_t));
  precond->context = context;
  precond->compute_diagonal = compute_diagonal;
  precond->dtor = dtor;
  precond->num_block_rows = num_block_rows;
  precond->block_size = block_size;
  precond->D = malloc(sizeof(real_t) * num_block_rows * block_size * block_size);
  memset(precond->D, 0, sizeof(real_t) * num_block_rows * block_size * block_size);

  preconditioner_vtable vtable = {.setup = block_jacobi_preconditioner_setup,
                                  .solve = block_jacobi_preconditioner_solve,
                                  .fprintf = block_jacobi_preconditioner_fprintf,
                                  .dtor = block_jacobi_preconditioner_dtor};
  return preconditioner_new("Block Jacobi preconditioner", precond, vtable);
}

void* block_jacobi_preconditioner_context(preconditioner_t* bj_precond)
{
  block_jacobi_preconditioner_t* precond = preconditioner_context(bj_precond);
  return precond->context;
}

