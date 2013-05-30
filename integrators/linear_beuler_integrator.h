#ifndef POLYMEC_LINEAR_BEULER_INTEGRATOR_H
#define POLYMEC_LINEAR_BEULER_INTEGRATOR_H

#include "integrators/integrator.h"
#include "core/sparse_lin_solver.h"

// Creates an integrator that uses the 1st-order implicit backward Euler 
// method to integrate a linear system of ordinary differential equations.
// The integrator uses a Krylov-powered linear solver under the hood.
// Arguments:
//   comm - The MPI communicator to use for parallel exchanges.
//   context - The context used to compute the properties of the linear 
//             operator A and the right-hand-side vector b, as well as 
//             any preconditioner.
//   compute_Ax - A function that applies the linear operator A to the solution vector x.
//   compute_rhs - A function that computes b at a given time.
//   apply_bcs - A function that applies discrete boundary conditions to the solution vector at a given time.
//   dtor - A destructor for freeing the context pointer (if any).
integrator_t* linear_beuler_integrator_new(MPI_Comm comm,
                                           void* context, 
                                           integrator_compute_Ax_func compute_Ax,
                                           integrator_compute_rhs_func compute_rhs,
                                           integrator_apply_bcs_func apply_bcs,
                                           integrator_dtor dtor);

// This function returns the sparse linear solver used by the backward Euler integrator.
sparse_lin_solver_t* linear_beuler_integrator_get_solver(integrator_t* integrator);

#endif

