// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/cpr_differencer.h"
#include "slu_ddefs.h"
#include "slu_util.h"
#include "core/linear_algebra.h"

struct cpr_differencer_t
{
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*F_dae)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* F_context;
};

cpr_differencer_t* cpr_differencer_new(int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                       int (*F_dae)(void* context, real_t, real_t* x, real_t* Fval),
                                       void* F_context,
                                       adj_graph_t* sparsity,
                                       int num_local_block_rows,
                                       int num_remote_block_rows,
                                       int block_size)
{
  cpr_differencer_t* diff = polymec_malloc(sizeof(cpr_differencer_t));
  return diff;
}

void cpr_differencer_free(cpr_differencer_t* diff)
{
  polymec_free(diff);
}

void cpr_differencer_compute(cpr_differencer_t* diff, 
                             real_t scale_factor,   
                             local_matrix_t* matrix);
{
}

struct local_matrix_t
{
  char* name;
  void* context;
  local_matrix_vtable vtable;
};

local_matrix_t* local_matrix_new(const char* name,
                                 void* context,
                                 local_matrix_vtable vtable)
{
  ASSERT(vtable->zero != NULL);
  ASSERT(vtable->add_column_vector != NULL);
  ASSERT(vtable->add_identity != NULL);
  ASSERT(vtable->solve != NULL);
  ASSERT(vtable->fprintf != NULL);
  local_matrix_t* matrix = polymec_malloc(sizeof(local_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  return matrix;
}
 
void local_matrix_free(local_matrix_t* matrix)
{
  if ((matrix->dtor != NULL) && (matrix->context != NULL))
    matrix->dtor(matrix->context);
  string_free(matrix->name);
  polymec_free(matrix);
}

char* local_matrix_name(local_matrix_t* matrix)
{
  return matrix->name;
}

void local_matrix_zero(local_matrix_t* matrix)
{
  matrix->vtable.zero(matrix->context);
}

void local_matrix_add_column_vector(local_matrix_t* matrix, 
                                    real_t scale_factor,
                                    int column,
                                    real_t* column_vector)
{
  matrix->vtable.add_column_vector(matrix->context,
                                   scale_factor,
                                   column,
                                   column_vector);
}

void local_matrix_add_identity(local_matrix_t* matrix, 
                               real_t scale_factor)
{
  matrix->vtable.add_identity(matrix->context, scale_factor);
}

bool local_matrix_solve(local_matrix_t* matrix, 
                       real_t* B,
                       real_t* x)
{
  return matrix->vtable.solve(matrix->context, B, x);
}

void local_matrix_fprintf(local_matrix_t* matrix, FILE* stream)
{
  matrix->vtable.fprintf(matrix->context, stream);
}

typedef struct 
{
  int num_block_rows;
  int *D_offsets, *B_offsets; // For variable block sizes.
  real_t* D;
  bool constant_block_size;
} bdm_t;

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
  A->constant_block_size = true;
  int bs0 = -1;
  for (int i = 0; i < num_block_rows; ++i)
  {
    int bs = block_sizes[i];
    if (bs0 == -1)
      bs0 = bs;
    else if (bs != bs0)
      A->constant_block_size = false;
    ASSERT(bs >= 1);
    A->D_offsets[i+1] = A->D_offsets[i] + bs*bs;
    A->B_offsets[i+1] = A->B_offsets[i] + bs;
  }
  int N = A->D_offsets[A->num_block_rows];
  A->D = polymec_malloc(sizeof(real_t) * N);

  char name[1024];
  if (A->constant_block_size)
    snprintf(name, 1024, "Block diagonal matrix (bs = %d)", bs0);
  else
    snprintf(name, 1024, "Variable block diagonal matrix");
  local_matrix_vtable vtable = {.dtor = bdm_dtor,
                                .zero = bdm_zero,
                                .add_identity = bdm_add_identity,
                                .add_column_vector = bdm_column_vector,
                                .solve = bdm_solve,
                                .fprintf = bdm_fprintf};
  return local_matrix_new(name, A, vtable);
}

local_matrix_t* sparse_local_matrix_new(adj_graph_t* sparsity)
{
}

