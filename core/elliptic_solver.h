#ifndef POLYMEC_ELLIPTIC_SOLVER_H
#define POLYMEC_ELLIPTIC_SOLVER_H

#include "core/polymec.h"
#include "core/table.h"
#include "core/index_space.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class solves the (time-independent) elliptic equation L(U) = S,
// where L is some linear operator, S is some source term, and U is the 
// the solution of the equation.
typedef struct elliptic_solver_t elliptic_solver_t;

// A function for computing the operator matrix at time t.
typedef void (*elliptic_solver_compute_op_mat_func)(void*, double_table_t*, double);

// A function for computing the source vector at time t. 
typedef void (*elliptic_solver_compute_source_func)(void*, double*, double);

// A function for applying boundary conditions to a linear system so 
// that the solution respects these boundary conditions. Arguments:
// context       - A pointer to the context used to apply boundary conditions.
// A             - The matrix in the linear system.
// b             - The right hand side in the linear system.
// t             - The time at which the boundary conditions are applied.
typedef void (*elliptic_solver_apply_bcs_func)(void*, double_table_t*, double*, double);

// A destructor function for the context object (if any).
typedef void (*elliptic_solver_dtor)(void*);

// This virtual table must be implemented by any elliptic_solver.
typedef struct 
{
  elliptic_solver_compute_op_mat_func    compute_operator_matrix;
  elliptic_solver_compute_source_func    compute_source_vector;
  elliptic_solver_apply_bcs_func         apply_bcs;
  elliptic_solver_dtor                   dtor;
} elliptic_solver_vtable;

// Creates a diffusion solver with the given name, context, and virtual table.
elliptic_solver_t* elliptic_solver_new(const char* name, 
                                       void* context,
                                       elliptic_solver_vtable vtable,
                                       index_space_t* index_space);

// Frees a diffusion solver.
void elliptic_solver_free(elliptic_solver_t* solver);

// Returns the name of the diffuson solver (internally stored).
char* elliptic_solver_name(elliptic_solver_t* solver);

// Returns the context object for this solver.
void* elliptic_solver_context(elliptic_solver_t* solver);

// Solves the elliptic equation at time t, placing the solution into 
// the given solution array.
void elliptic_solver_solve(elliptic_solver_t* solver,
                           double t,
                           double* solution);

#ifdef __cplusplus
}
#endif

#endif

