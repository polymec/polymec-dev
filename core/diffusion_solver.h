#ifndef POLYMEC_DIFFUSION_SOLVER_H
#define POLYMEC_DIFFUSION_SOLVER_H

#include "core/polymec.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class solves the linear diffusion equation dU/dt = D(U) + S, 
// where U is some quantity undergoing diffusion, D is a diffusion operator, 
// and S is a source.
typedef struct diffusion_solver_t diffusion_solver_t;

// A function for creating a matrix.
typedef void (*diffusion_solver_create_mat_func)(void*, HYPRE_IJMatrix*);

// A function for creating a vector.
typedef void (*diffusion_solver_create_vec_func)(void*, HYPRE_IJVector*);

// A function for creating a linear solver.
typedef void (*diffusion_solver_create_ksp_func)(void*, HYPRE_Solver*);

// A function for computing the diffusion matrix at time t.
typedef void (*diffusion_solver_compute_diff_mat_func)(void*, HYPRE_IJMatrix, double);

// A function for computing the source vector at time t. 
typedef void (*diffusion_solver_compute_source_func)(void*, HYPRE_IJVector, double);

// A function for applying boundary conditions to a linear system so 
// that the solution respects these boundary conditions. Arguments:
// void* context - A pointer to the context used to apply boundary conditions.
// HYPRE_IJMatrix A         - The matrix in the linear system.
// HYPRE_IJVector b         - The right hand side in the linear system.
// t             - The time at which the boundary conditions are applied.
typedef void (*diffusion_solver_apply_bcs_func)(void*, HYPRE_IJMatrix, HYPRE_IJVector, double);

// A destructor function for the context object (if any).
typedef void (*diffusion_solver_dtor)(void*);

// This virtual table must be implemented by any diffusion_solver.
typedef struct 
{
  diffusion_solver_create_mat_func        create_matrix;
  diffusion_solver_create_vec_func        create_vector;
  diffusion_solver_create_ksp_func        create_ksp;
  diffusion_solver_compute_diff_mat_func  compute_diffusion_matrix;
  diffusion_solver_compute_source_func    compute_source_vector;
  diffusion_solver_apply_bcs_func         apply_bcs;
  diffusion_solver_dtor                   dtor;
} diffusion_solver_vtable;

// Creates a diffusion solver with the given name, context, and virtual table.
diffusion_solver_t* diffusion_solver_new(const char* name, 
                                         void* context,
                                         diffusion_solver_vtable vtable);

// Frees a diffusion solver.
void diffusion_solver_free(diffusion_solver_t* solver);

// Returns the name of the diffuson solver (internally stored).
char* diffusion_solver_name(diffusion_solver_t* solver);

// Returns the context object for this solver.
void* diffusion_solver_context(diffusion_solver_t* solver);

// Solves the diffusion equation using the L-stable, first-order-accurate 
// backward Euler method. Arguments:
// solver - the diffusion solver.
// t1           - The start time for the integration.
// sol1         - An array containing the solution at time t1.
// t2           - The end time for the integration.
// sol2         - An array that will store the solution at time t2.
void diffusion_solver_euler(diffusion_solver_t* solver,
                            double t1, double* sol1, 
                            double t2, double* sol2);

// This function solves the diffusion equation using the "TGA" algorithm 
// presented by Twizell, et al (1996). Arguments:
// t1           - The start time for the integration.
// sol1         - An array containing the solution at time t1.
// t2           - The end time for the integration.
// sol2         - An array that will store the solution at time t2.
void diffusion_solver_tga(diffusion_solver_t* solver,
                          double t1, double* sol1,
                          double t2, double* sol2);

#ifdef __cplusplus
}
#endif

#endif

