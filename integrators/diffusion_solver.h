#ifndef POLYMEC_DIFFUSION_SOLVER_H
#define POLYMEC_DIFFUSION_SOLVER_H

#include "core/polymec.h"
#include "core/table.h"
#include "core/index_space.h"

// This class solves the linear diffusion equation dU/dt = D(U) + S, 
// where U is some quantity undergoing diffusion, D is a diffusion operator, 
// and S is a source.
typedef struct diffusion_solver_t diffusion_solver_t;

// A function for computing the coefficients of the diffusion matrix at time t.
typedef void (*diffusion_solver_compute_diff_mat_func)(void*, double_table_t*, double);

// A function for computing the components of the source vector at time t. 
typedef void (*diffusion_solver_compute_source_func)(void*, double*, double);

// A function for applying boundary conditions to a linear system so 
// that the solution respects these boundary conditions. Arguments:
// void* context - A pointer to the context used to apply boundary conditions.
// A             - The coefficients of the matrix in the linear system.
// b             - The components of the right hand side in the linear system.
// t             - The time at which the boundary conditions are applied.
typedef void (*diffusion_solver_apply_bcs_func)(void*, double_table_t*, double*, double);

// A destructor function for the context object (if any).
typedef void (*diffusion_solver_dtor)(void*);

// This virtual table must be implemented by any diffusion_solver.
typedef struct 
{
  diffusion_solver_compute_diff_mat_func  compute_diffusion_matrix;
  diffusion_solver_compute_source_func    compute_source_vector;
  diffusion_solver_apply_bcs_func         apply_bcs;
  diffusion_solver_dtor                   dtor;
} diffusion_solver_vtable;

// Creates a diffusion solver with the given name, context, and virtual table.
// Also needs an index space for constructing the linear system.
diffusion_solver_t* diffusion_solver_new(const char* name, 
                                         void* context,
                                         diffusion_solver_vtable vtable,
                                         index_space_t* index_space);

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

// This function solves the diffusion equation using the A-stable "TGA" 
// algorithm presented by Twizell, et al (1996). Arguments:
// t1           - The start time for the integration.
// sol1         - An array containing the solution at time t1.
// t2           - The end time for the integration.
// sol2         - An array that will store the solution at time t2.
void diffusion_solver_tga(diffusion_solver_t* solver,
                          double t1, double* sol1,
                          double t2, double* sol2);

// This function solves the diffusion equation using the Crank-Nicolson
// algorithm, which is 2nd-order, but not A-stable.
// t1           - The start time for the integration.
// sol1         - An array containing the solution at time t1.
// t2           - The end time for the integration.
// sol2         - An array that will store the solution at time t2.
void diffusion_solver_crank_nicolson(diffusion_solver_t* solver,
                                     double t1, double* sol1,
                                     double t2, double* sol2);

#endif

