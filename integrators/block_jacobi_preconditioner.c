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

#include "integrators/block_jacobi_preconditioner.h"

typedef struct 
{
  int (*F)(void* context, real_t t, real_t* x, real_t* F);
  void (*communicate)(void* context, real_t t, real_t* x);
  void* context;

  int num_block_rows;
  int block_size; 

} block_jacobi_preconditioner_t;

// Transforms A to (I + gamma * A).
static void bd_scale_and_shift(void* context, real_t gamma)
{
}

// Returns the (i, j)th entry in a block-diagonal.
static real_t bd_coeff(void* context, int i, int j)
{
  // FIXME
  return 0.0;
}

static void bd_dtor(void* context)
{
  // FIXME
}

static preconditioner_matrix_t* block_jacobi_preconditioner_matrix(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  preconditioner_matrix_vtable vtable = {.scale_and_shift = bd_scale_and_shift,
                                         .coeff = bd_coeff,
                                         .dtor = bd_dtor};
  int num_rows = precond->num_block_rows * precond->block_size;
  return preconditioner_matrix_new("Block-diagonal", NULL, vtable, num_rows);
}

static void block_jacobi_preconditioner_compute_jacobian(void* context, real_t t, real_t* x, preconditioner_matrix_t* mat)
{
  block_jacobi_preconditioner_t* precond = context;
}

static bool block_jacobi_preconditioner_solve(void* context, preconditioner_matrix_t* A, real_t* B)
{
  block_jacobi_preconditioner_t* precond = context;
  // FIXME
  return true;
}

static void block_jacobi_preconditioner_dtor(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
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

  preconditioner_vtable vtable = {.matrix = block_jacobi_preconditioner_matrix,
                                  .compute_jacobian = block_jacobi_preconditioner_compute_jacobian,
                                  .solve = block_jacobi_preconditioner_solve,
                                  .dtor = block_jacobi_preconditioner_dtor};
  return preconditioner_new("Block Jacobi preconditioner", precond, vtable);
}
