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

#include "integrators/preconditioner.h"

struct preconditioner_t 
{
  char* name;
  void* context;
  preconditioner_vtable vtable;
};

struct preconditioner_matrix_t 
{
  char* name;
  void* context;
  preconditioner_matrix_vtable vtable;
};

preconditioner_t* preconditioner_new(const char* name,
                                     void* context,
                                     preconditioner_vtable vtable)
{
  preconditioner_t* precond = malloc(sizeof(preconditioner_t));
  precond->name = string_dup(name);
  precond->context = context;
  precond->vtable = vtable;

  return precond;
}

void preconditioner_free(preconditioner_t* precond)
{
  if ((precond->vtable.dtor != NULL) && (precond->context != NULL))
    precond->vtable.dtor(precond->context);
  free(precond->name);
  free(precond);
}

char* preconditioner_name(preconditioner_t* precond)
{
  return precond->name;
}

void* preconditioner_context(preconditioner_t* precond)
{
  return precond->context;
}

preconditioner_matrix_t* preconditioner_matrix(preconditioner_t* precond)
{
  return precond->vtable.matrix(precond->context);
}

void preconditioner_compute_jacobian(preconditioner_t* precond,
                                     real_t t,
                                     real_t* x,
                                     preconditioner_matrix_t* mat)
{
  precond->vtable.compute_jacobian(precond->context, t, x, mat);
}

void preconditioner_solve(preconditioner_t* precond,
                          preconditioner_matrix_t* mat,
                          real_t* rhs)
{
  precond->vtable.solve(precond->context, mat, rhs);
}

preconditioner_matrix_t* preconditioner_matrix_new(const char* name,
                                                   void* context,
                                                   preconditioner_matrix_vtable vtable,
                                                   int num_rows)
{
  ASSERT(vtable.scale_and_shift != NULL);
  ASSERT(vtable.coeff != NULL);
  preconditioner_matrix_t* mat = malloc(sizeof(preconditioner_matrix_t));
  mat->name = string_dup(name);
  mat->context = context;
  mat->vtable = vtable;

  return mat;
}

void preconditioner_matrix_free(preconditioner_matrix_t* mat)
{
  if ((mat->vtable.dtor != NULL) && (mat->context != NULL))
    mat->vtable.dtor(mat->context);
  free(mat->name);
  free(mat);
}

void* preconditioner_matrix_context(preconditioner_matrix_t* mat)
{
  return mat->context;
}

void preconditioner_matrix_scale_and_shift(preconditioner_matrix_t* mat, real_t gamma)
{
  mat->vtable.scale_and_shift(mat->context, gamma);
}

real_t preconditioner_matrix_coeff(preconditioner_matrix_t* mat, int i, int j)
{
  return mat->vtable.coeff(mat->context, i, j);
}

