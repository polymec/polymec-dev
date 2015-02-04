// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/local_matrix.h"

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
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_column_vector != NULL);
  ASSERT(vtable.add_identity != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(vtable.fprintf != NULL);
  local_matrix_t* matrix = polymec_malloc(sizeof(local_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  return matrix;
}
 
void local_matrix_free(local_matrix_t* matrix)
{
  if ((matrix->vtable.dtor != NULL) && (matrix->context != NULL))
    matrix->vtable.dtor(matrix->context);
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

