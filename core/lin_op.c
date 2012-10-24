#include <stdlib.h>
#include <string.h>
#include <gc/gc.h>
#include "core/lin_op.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  free(op->name);
  free(op);
}

lin_op_t* lin_op_new(const char* name, void* context, lin_op_vtable vtable,
                     mesh_t* mesh)
{
  ASSERT(vtable.stencil_size != NULL);
  ASSERT(vtable.compute_offsets != NULL);
  ASSERT(vtable.compute_weights != NULL);
  lin_op_t* L = GC_MALLOC(sizeof(lin_op_t));
  L->name = strdup(name);
  L->context = context;
  L->vtable = vtable;
  L->mesh = mesh;
  GC_register_finalizer(L, &lin_op_free, L, NULL, NULL);
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

void lin_op_compute_offsets(lin_op_t* op, int index, int* offsets)
{
  ASSERT(index >= 0);
  ASSERT(offsets != NULL);
  return op->vtable.compute_offsets(op->context, op->mesh, index, offsets);
}

void lin_op_compute_weights(lin_op_t* op, int index, int* offsets, double* weights)
{
  ASSERT(index >= 0);
  ASSERT(offsets != NULL);
  ASSERT(weights != NULL);
  return op->vtable.compute_weights(op->context, op->mesh, index, offsets, weights);
}

#ifdef __cplusplus
}
#endif

