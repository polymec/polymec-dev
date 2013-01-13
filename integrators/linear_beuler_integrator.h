#ifndef POLYMEC_LINEAR_BEULER_INTEGRATOR_H
#define POLYMEC_LINEAR_BEULER_INTEGRATOR_H

#include "core/integrator.h"
#include "sundials/sundials_iterative.h"

#ifdef __cplusplus
extern "C" {
#endif

// A prototype for a function that computes right-hand-side vectors for 
// linear backward Euler integrators at a given time.
typedef void (*linear_beuler_compute_rhs)(void*, double, N_Vector);

// Creates an integrator that uses the 1st-order implicit backward Euler 
// method to integrate a linear system of ordinary differential equations.
// This version solves the linear system using the Generalized Minimum 
// Residual (GMRES) method. Arguments:
//   context - The context used to compute the properties of the linear 
//             operator A and the right-hand-side vector b, as well as 
//             any preconditioner.
//   Ax - A function that applies A to the vector x.
//   compute_rhs - A function that computes b at a given time.
//   precond_type - A flag that determines the type of preconditioning to 
//                  be used. Possible values are PREC_NONE for no preconditioning,
//                  PREC_LEFT for left preconditioning, PREC_RIGHT for 
//                  right preconditioning, or PREC_BOTH for both left and right.
//   precond - A function that solves the preconditioner system Pz = r. 
//             This can be NULL if precond_type is PREC_NONE.
//   gram_schmidt - A flag that determines whether to use the classical 
//                  (CLASSICAL_GS) or the modified (MODIFIED_GS) Gram-Schmidt 
//                  orthogonalization process.
//   delta - The tolerance on the L2 norm of the scaled, preconditioned 
//           residual.
//   max_restarts - The maximum number of times the algorithm is allowed 
//                  to restart.
//   dtor - A destructor for freeing the context pointer (if any).
integrator_t* gmres_linear_beuler_integrator_new(void* context, 
                                                 ATimesFn Ax,
                                                 linear_beuler_compute_rhs compute_rhs,
                                                 int precond_type,
                                                 PSolveFn precond,
                                                 int gram_schmidt,
                                                 double delta,
                                                 int max_restarts,
                                                 integrator_dtor dtor);

// Creates an integrator that uses the 1st-order implicit backward Euler 
// method to integrate a linear system of ordinary differential equations.
// This version solves the linear system using the stabilized Bi-conjugate
// Gradient (Bi-CGSTAB) method. Arguments:
//   context - The context used to compute the properties of the linear 
//             operator A and the right-hand-side vector b, as well as 
//             any preconditioner.
//   Ax - A function that applies A to the vector x.
//   compute_rhs - A function that computes b at a given time.
//   precond_type - A flag that determines the type of preconditioning to 
//                  be used. Possible values are PREC_NONE for no preconditioning,
//                  PREC_LEFT for left preconditioning, PREC_RIGHT for 
//                  right preconditioning, or PREC_BOTH for both left and right.
//   precond - A function that solves the preconditioner system Pz = r. 
//             This can be NULL if precond_type is PREC_NONE.
//   delta - The tolerance on the L2 norm of the scaled, preconditioned 
//           residual.
//   dtor - A destructor for freeing the context pointer (if any).
integrator_t* bicgstab_linear_beuler_integrator_new(void* context, 
                                                    ATimesFn Ax,
                                                    linear_beuler_compute_rhs compute_rhs,
                                                    int precond_type,
                                                    PSolveFn precond,
                                                    double delta,
                                                    integrator_dtor dtor);

// This function retrieves information pertaining to the most recent step
// performed by the given linear backward Euler integrator.
void linear_beuler_integrator_get_info(integrator_t* integrator,
                                       double* res_l2_norm,
                                       int* num_linear_iterations,
                                       int* num_precond_solves);

#ifdef __cplusplus
}
#endif

#endif

