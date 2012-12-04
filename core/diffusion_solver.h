#ifndef POLYMEC_SOLVE_DIFFUSION_EQ_H
#define POLYMEC_SOLVE_DIFFUSION_EQ_H

#include "core/polymec.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class solves the linear diffusion equation dU/dt = D(U) + S, 
// where U is some quantity undergoing diffusion, D is a diffusion operator, 
// and S is a source.
typedef struct diffusion_solver_t diffusion_solver_t;

// Functions of this type apply boundary conditions to a linear system so 
// that the solution respects these boundary conditions. Arguments:
// void* context - A pointer to the context used to apply boundary conditions.
// Mat A         - The matrix in the linear system.
// Vec b         - The right hand side in the linear system.
// t             - The time at which the boundary conditions are applied.
typedef void (*diffusion_solver_apply_bcs_func)(void*, Mat, Vec, double);

// Creates a diffusion solver on the given MPI communicator, with the 
// given function that is used to apply diffusion boundary conditions to 
// a diffusion linear system Ax = b (where x is the solution. Also provided 
// are a context pointer passed to apply_bcs, and a destructor for that pointer.
diffusion_solver_t* diffusion_solver_new(MPI_Comm comm,
                                         diffusion_solver_apply_bcs_func apply_bcs, 
                                         void* context,
                                         void (*context_dtor)(void*));

// Frees a diffusion solver.
void diffusion_solver_free(diffusion_solver_t* solver);

// Solves the diffusion equation using the L-stable, first-order-accurate 
// backward Euler method. Arguments:
// solver - the diffusion solver.
// diffusion_op - A matrix representing the diffusion operator D.
// source       - A vector representing the source term S.
// t1           - The start time for the integration.
// t2           - The end time for the integration.
// sol1         - A vector containing the solution at time t1.
// sol2         - A vector that will store the solution at time t2.
void diffusion_solver_euler(diffusion_solver_t* solver,
                            Mat diffusion_op, 
                            Vec source, 
                            double t1,
                            double t2,
                            Vec sol1, 
                            Vec sol2);

// This function solves the diffusion equation using the "TGA" algorithm 
// presented by Twizell, et al (1996). Arguments:
// diffusion_op - A matrix representing the diffusion operator D.
// source1      - A vector representing the source term S at the start time.
// source2      - A vector representing the source term S at the end time.
// t1           - The start time for the integration.
// t2           - The end time for the integration.
// sol1         - A vector containing the solution at time t1.
// sol2         - A vector that will store the solution at time t2.
void diffusion_solver_tga(diffusion_solver_t* solver,
                          Mat diffusion_op, 
                          Vec source1, 
                          Vec source2, 
                          double t1,
                          double t2,
                          Vec sol1, 
                          Vec sol2);

#ifdef __cplusplus
}
#endif

#endif

