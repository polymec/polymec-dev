#ifndef POLYMEC_SPARSE_LIN_SOLVER_H
#define POLYMEC_SPARSE_LIN_SOLVER_H

#include "sundials/sundials_iterative.h"
#include "core/polymec.h"
#if USE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

// These macros help with manipulating serial and parallel N_Vector objects.
#if USE_MPI
#define NV_DATA(v) NV_DATA_P(v)
#define NV_LOCLENGTH(v) NV_LOCLENGTH_P(v)
#define NV_GLOBLENGTH(v) NV_GLOBLENGTH_P(v)
#define NV_Ith(v) NV_Ith_P(v)
#else
#define NV_DATA(v) NV_DATA_S(v)
#define NV_LOCLENGTH(v) NV_LENGTH_S(v)
#define NV_GLOBLENGTH(v) NV_LENGTH_S(v)
#define NV_Ith(v) NV_Ith_S(v)
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
// subclass. Only solve is strictly required.
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

