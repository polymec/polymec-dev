// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string.h>
#include <gc/gc.h>
#include "core/lin_op.h"

struct lin_op_t
{
  char* name;
  void* context;
  lin_op_vtable vtable;
  mesh_t* mesh;
};

static void lin_op_free(void* ctx, void* dummy)
{
  lin_op_t* op = (lin_op_t*)ctx;
  if (op->vtable.dtor)
    op->vtable.dtor(op->context);
  op->context = NULL;
  free(op->name);
  free(op);
}

lin_op_t* lin_op_new(const char* name, void* context, lin_op_vtable vtable,
                     mesh_t* mesh)
{
  ASSERT(vtable.stencil_size != NULL);
  ASSERT(vtable.compute_stencil != NULL);
  lin_op_t* L = GC_MALLOC(sizeof(lin_op_t));
  L->name = strdup(name);
  L->context = context;
  L->vtable = vtable;
  L->mesh = mesh;
  GC_register_finalizer(L, lin_op_free, L, NULL, NULL);
  return L;
}

char* lin_op_name(lin_op_t* op)
{
  return op->name;
}

void* lin_op_context(lin_op_t* op)
{
  return op->context;
}

int lin_op_stencil_size(lin_op_t* op, int index)
{
  ASSERT(index >= 0);
  return op->vtable.stencil_size(op->context, op->mesh, index);
}

void lin_op_compute_stencil(lin_op_t* op, int index, int* offsets, double* weights)
{
  ASSERT(index >= 0);
  ASSERT(offsets != NULL);
  ASSERT(weights != NULL);
  op->vtable.compute_stencil(op->context, op->mesh, index, offsets, weights);
}

void lin_op_apply(lin_op_t* op, double* field, double* Lfield)
{
  ASSERT(field != NULL);
  ASSERT(Lfield != NULL);
  op->vtable.apply(op->context, op->mesh, field, Lfield);
}

