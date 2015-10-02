// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/local_matrix.h"
#include "core/timer.h"

struct local_matrix_t
{
  char* name;
  void* context;
  local_matrix_vtable vtable;
  int num_rows;
};

local_matrix_t* local_matrix_new(const char* name,
                                 void* context,
                                 local_matrix_vtable vtable,
                                 int num_rows)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_column_vector != NULL);
  ASSERT(vtable.add_row_vector != NULL);
  ASSERT(vtable.add_identity != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(vtable.fprintf != NULL);
  ASSERT(vtable.value != NULL);
  ASSERT(vtable.set_value != NULL);
  ASSERT(vtable.get_diag != NULL);
  ASSERT(vtable.matvec != NULL);
  ASSERT(num_rows > 0);
  local_matrix_t* matrix = polymec_malloc(sizeof(local_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  matrix->num_rows = num_rows;
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

void* local_matrix_context(local_matrix_t* matrix)
{
  return matrix->context;
}

int local_matrix_num_rows(local_matrix_t* matrix)
{
  return matrix->num_rows;
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
  START_FUNCTION_TIMER();
  matrix->vtable.add_column_vector(matrix->context,
                                   scale_factor,
                                   column,
                                   column_vector);
  STOP_FUNCTION_TIMER();
}

void local_matrix_add_row_vector(local_matrix_t* matrix, 
                                 real_t scale_factor,
                                 int row,
                                 real_t* row_vector)
{
  START_FUNCTION_TIMER();
  matrix->vtable.add_row_vector(matrix->context,
                                scale_factor,
                                row,
                                row_vector);
  STOP_FUNCTION_TIMER();
}

void local_matrix_add_identity(local_matrix_t* matrix, 
                               real_t scale_factor)
{
  START_FUNCTION_TIMER();
  matrix->vtable.add_identity(matrix->context, scale_factor);
  STOP_FUNCTION_TIMER();
}

bool local_matrix_solve(local_matrix_t* matrix, real_t* B, real_t* x)
{
  START_FUNCTION_TIMER();
  bool result = matrix->vtable.solve(matrix->context, B, x);
  STOP_FUNCTION_TIMER();
  return result;
}

void local_matrix_fprintf(local_matrix_t* matrix, FILE* stream)
{
  matrix->vtable.fprintf(matrix->context, stream);
}

real_t local_matrix_value(local_matrix_t* matrix, int i, int j)
{
  ASSERT(i >= 0);
  ASSERT(j >= 0);
  return matrix->vtable.value(matrix->context, i, j);
}

void local_matrix_set_value(local_matrix_t* matrix, int i, int j, real_t value)
{
  ASSERT(i >= 0);
  ASSERT(j >= 0);
  matrix->vtable.set_value(matrix->context, i, j, value);
}

void local_matrix_get_diagonal(local_matrix_t* matrix, real_t* diag)
{
  START_FUNCTION_TIMER();
  matrix->vtable.get_diag(matrix->context, diag);
  STOP_FUNCTION_TIMER();
}

void local_matrix_matvec(local_matrix_t* matrix, real_t* x, real_t* Ax)
{
  START_FUNCTION_TIMER();
  matrix->vtable.matvec(matrix->context, x, Ax);
  STOP_FUNCTION_TIMER();
}
