#ifndef POLYMEC_SPARSE_LIN_SOLVER_H
#define POLYMEC_SPARSE_LIN_SOLVER_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

// This type represents the outcome of an iterative sparse linear solve.
typedef enum
{
  SPARSE_LIN_SOLVER_CONVERGED,
  SPARSE_LIN_SOLVER_REDUCED_RESIDUAL,
  SPARSE_LIN_SOLVER_DID_NOT_CONVERGE,
  SPARSE_LIN_SOLVER_FAILED_RECOVERABLY,
  SPARSE_LIN_SOLVER_FAILED_UNRECOVERABLY
} sparse_lin_solver_outcome_t;

// This is the maximum length of the string that contains any diagnostic details about the 
// outcome of the most recent linear solve.
#define SPARSE_LIN_SOLVER_OUTCOME_DETAILS_MAXLEN 2048

// This class provides an interface for a sparse iterative integrator that 
// solves a system of linear equations using a scaled, preconditioned Krylov 
// subsequence method.
typedef struct sparse_lin_solver_t sparse_lin_solver_t;

// A prototype for a function that solves a linear system.
typedef sparse_lin_solver_outcome_t (*sparse_lin_solver_solve_func)(void* context, void* x, void* b, double* res_norm, int* num_lin_iterations, int* num_precond_solves, char* outcome_details);

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

// Solves the linear system of equations for the given right hand side vector b,
// placing the solution in the vector x. Note that the vectors are represented
// in a "neutral" void* form.
sparse_lin_solver_outcome_t sparse_lin_solver_solve(sparse_lin_solver_t* solver, void* x, void* b);

// This function retrieves information pertaining to the most recent 
// successful solve performed by the given solver.
void sparse_lin_solver_get_info(sparse_lin_solver_t* solver,
                                double* res_l2_norm,
                                int* num_linear_iterations,
                                int* num_precond_solves);

// This retrieves an internally-stored string that contains diagnostic details about the 
// outcome (usually unsuccessful) of the most recent linear solve.
char* sparse_lin_solver_outcome_details(sparse_lin_solver_t* solver);

#ifdef __cplusplus
}
#endif

#endif

