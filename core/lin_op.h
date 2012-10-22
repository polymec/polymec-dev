#ifndef ARBI_LIN_OP_H
#define ARBI_LIN_OP_H

#include "core/arbi.h"
#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// A linear operator. Objects of this type are garbage-collected.
typedef struct lin_op_t lin_op_t;

// A function for determining the number of weights in the stencil of 
// the linear operator for a given cell.
typedef int (*lin_op_stencil_size_func)(void*, mesh_t*, cell_t*);

// A function for computing the offsets (in index space) for the weights
// in the stencil of the linear operator for a given cell.
typedef void (*lin_op_compute_offsets_func)(void*, mesh_t*, cell_t*, int*);

// A function for computing the offsets (in index space) for the weights
// in the stencil of the linear operator for a given cell.
typedef void (*lin_op_compute_weights_func)(void*, mesh_t*, cell_t*, int*, double*);

// A destructor function for the context object (if any).
typedef void (*lin_op_dtor)(void*);

// This virtual table must be implemented by any linear operator.
typedef struct 
{
  lin_op_stencil_size_func    stencil_size;
  lin_op_compute_offsets_func compute_offsets;
  lin_op_compute_weights_func compute_weights;
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
// applied to the given cell.
int lin_op_stencil_size(lin_op_t* op, cell_t* cell);

// Computes the offsets (in index space) for the weights in the stencil 
// of the linear operator applied to the given cell.
void lin_op_compute_offsets(lin_op_t* op, cell_t* cell, int* offsets);

// Computes the weights for the stencil of the linear operator applied to 
// the given cell.
void lin_op_compute_weights(lin_op_t* op, cell_t* cell, int* offsets, double* weights);

// These methods can be used in the virtual table of any 2nd-order finite 
// volume linear operator.
int fv2_lin_op_stencil_size(void* context, mesh_t* mesh, cell_t* cell);
void fv2_lin_op_compute_offsets(void* context, mesh_t* mesh, cell_t* cell, int* offsets);


#ifdef __cplusplus
}
#endif

#endif

