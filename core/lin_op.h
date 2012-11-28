#ifndef POLYMEC_LIN_OP_H
#define POLYMEC_LIN_OP_H

#include "core/polymec.h"
#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// A linear operator. Objects of this type are garbage-collected.
typedef struct lin_op_t lin_op_t;

// A function for determining the number of weights in the stencil of 
// the linear operator for a given cell.
typedef int (*lin_op_stencil_size_func)(void*, mesh_t*, int);

// A function for computing the offsets (in index space) and the weights
// for the stencil of the linear operator for a given cell.
typedef void (*lin_op_compute_stencil_func)(void*, mesh_t*, int, int*, double*);

// A function for directly applying the linear operator to a field.
typedef void (*lin_op_apply_func)(void*, mesh_t*, double*, double*);

// A destructor function for the context object (if any).
typedef void (*lin_op_dtor)(void*);

// This virtual table must be implemented by any linear operator.
typedef struct 
{
  lin_op_stencil_size_func    stencil_size;
  lin_op_compute_stencil_func compute_stencil;
  lin_op_apply_func           apply;
  lin_op_dtor                 dtor;
} lin_op_vtable;

// Creates an instance of a linear operator with the given name and vtable, 
// associating it with the given mesh.
lin_op_t* lin_op_new(const char* name, void* context, lin_op_vtable vtable, mesh_t* mesh);

// Returns the name of the linear operator.
char* lin_op_name(lin_op_t* op);

// Returns the context object associated with the model (if any).
void* lin_op_context(lin_op_t* op);

// Returns the number of weights in the stencil for the linear operator 
// applied to the geometric element with the given index.
int lin_op_stencil_size(lin_op_t* op, int index);

// Computes the offsets (in index space) and weights for the stencil of 
// the linear operator applied to the geometric element with the given index.
void lin_op_compute_stencil(lin_op_t* op, int index, int* offsets, double* weights);

// Applies the linear operator to the given field.
void lin_op_apply(lin_op_t* op, double* field, double* Lfield);

#ifdef __cplusplus
}
#endif

#endif

