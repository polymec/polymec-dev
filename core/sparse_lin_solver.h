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

// These macros and functions help with creating and manipulating serial 
// and parallel N_Vector objects.
#if USE_MPI
#define NV_DATA(v) NV_DATA_P(v)
#define NV_LOCLENGTH(v) NV_LOCLENGTH_P(v)
#define NV_GLOBLENGTH(v) NV_GLOBLENGTH_P(v)
#define NV_Ith(v) NV_Ith_P(v)

static inline N_Vector N_VNew(MPI_Comm comm, long int local_length)
{
  int N;
  MPI_Allreduce(&local_length, &N, 1, MPI_LONG, MPI_SUM, comm);
  return N_VNew_Parallel(comm, local_length, N);
}

static inline N_Vector N_VMake(MPI_Comm comm, long int local_length, double* v_data)
{
  int N;
  MPI_Allreduce(&local_length, &N, 1, MPI_LONG, MPI_SUM, comm);
  return N_VMake_Parallel(comm, local_length, N, v_data);
}

#else
#define NV_DATA(v) NV_DATA_S(v)
#define NV_LOCLENGTH(v) NV_LENGTH_S(v)
#define NV_GLOBLENGTH(v) NV_LENGTH_S(v)
#define NV_Ith(v) NV_Ith_S(v)

static inline N_Vector N_VNew(MPI_Comm comm, long int local_length)
{
  return N_VNew_Serial(local_length);
}

static inline N_Vector N_VMake(MPI_Comm comm, long int local_length, double* v_data)
{
  return N_VMake_Serial(local_length, v_data);
}

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
typedef sparse_lin_solver_outcome_t (*sparse_lin_solver_solve_func)(void* context, N_Vector x, N_Vector b, double* res_norm, int* num_lin_iterations, int* num_precond_solves, char* outcome_details);

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
sparse_lin_solver_outcome_t sparse_lin_solver_solve(sparse_lin_solver_t* solver, N_Vector x, N_Vector b);

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

