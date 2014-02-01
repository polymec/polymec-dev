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

#include "integrators/jacobi_preconditioners.h"

typedef struct
{
  void* context;
  int (*residual_func)(void* context, real_t t, real_t* x, real_t* F);
  void (*dtor)(void* context);
  int block_size;
} jacobi_precond_t;

static void jacobi_dtor(void* context)
{
  jacobi_precond_t* precond = context;
  if ((precond->context != NULL) && (precond->dtor != NULL))
    precond->dtor(precond->context);
}

static preconditioner_matrix_t* jacobi_matrix(void* context)
{
  // Bogus!
  return NULL;
}

static void jacobi_compute_jacobian(void* context, real_t t, real_t* x, preconditioner_matrix_t* mat)
{
  // Bogus!
}

static void jacobi_solve(void* context, preconditioner_matrix_t* A, real_t* B)
{
  // Bogus!
}

preconditioner_t* jacobi_preconditioner_new(void* context,
                                            int (*residual_func)(void* context, real_t t, real_t* x, real_t* F),
                                            void (*context_dtor),
                                            int block_size)
{
  ASSERT(block_size > 0);
  jacobi_precond_t* jacobi = malloc(sizeof(jacobi_precond_t));
  jacobi->context = context;
  jacobi->residual_func = residual_func;
  jacobi->block_size = block_size;

  preconditioner_vtable vtable = {.matrix = jacobi_matrix,
                                  .compute_jacobian = jacobi_compute_jacobian,
                                  .solve = jacobi_solve,
                                  .dtor = jacobi_dtor};
  return preconditioner_new("Jacobi preconditioner", jacobi, vtable);
}

void jacobi_preconditioner_get_data(preconditioner_t* precond,
                                    void** context,
                                    int (**residual_func)(void* context, real_t t, real_t* x, real_t* F),
                                    int* block_size)
{
  ASSERT(strcmp(preconditioner_name(precond), "Jacobi preconditioner") == 0);
  jacobi_precond_t* jacobi = preconditioner_context(precond);
  *context = jacobi->context;
  *residual_func = jacobi->residual_func;
  *block_size = jacobi->block_size;
}

