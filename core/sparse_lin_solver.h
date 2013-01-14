#ifndef POLYMEC_SPARSE_LIN_SOLVER_H
#define POLYMEC_SPARSE_LIN_SOLVER_H

#include "sundials/sundials_iterative.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class provides an interface for a sparse iterative integrator that 
// solves a system of linear equations using a scaled, preconditioned Krylov 
// subsequence method.
typedef struct sparse_lin_solver_t sparse_lin_solver_t;

// A prototype for a function that solves a linear system.
typedef void (*sparse_lin_solver_solve_func)(void* context, N_Vector x, N_Vector b, double* res_norm, int* num_lin_iterations, int* num_precond_solves);

// A destructor function for the context object (if any).
typedef void (*sparse_lin_solver_dtor)(void*);

// This virtual table must be implemented for each sparse linear solver 
// subclass.
typedef struct
{
  sparse_lin_solver_solve_func solve;
  sparse_lin_solver_dtor dtor;
} sparse_lin_solver_vtable;

// Creates new sparse linear solver with the given name, context, and virtual 
// table.
sparse_lin_solver_t* sparse_lin_solver_new(const char* name, void* context,
                                           sparse_lin_solver_vtable vtable);

// Frees the storage associated with the given linear solver.
void sparse_lin_solver_free(sparse_lin_solver_t* solver);

// Solves the linear system of equations for the given right hand side vector.
void sparse_lin_solver_solve(sparse_lin_solver_t* solver, N_Vector x, N_Vector b);

// This function retrieves information pertaining to the most recent 
// solve performed by the given solver.
void sparse_lin_solver_get_info(sparse_lin_solver_t* solver,
                                double* res_l2_norm,
                                int* num_linear_iterations,
                                int* num_precond_solves);

#ifdef __cplusplus
}
#endif

#endif

