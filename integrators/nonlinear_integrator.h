#ifndef POLYMEC_NONLINEAR_INTEGRATOR_H
#define POLYMEC_NONLINEAR_INTEGRATOR_H

#include "core/integrator.h"
#include "cvode/cvode.h"
#include "cvode/cvode_spils.h"

#ifdef __cplusplus
extern "C" {
#endif

// Preconditioner setup function.
typedef CVSpilsPrecSetupFn nonlinear_integrator_precond_setup_func;

// Preconditioner solve function.
typedef CVSpilsPrecSolveFn nonlinear_integrator_precond_solve_func;

// Creates a integrator that uses CVODE to integrate a nonlinear system 
// of equations.
//   comm - The MPI communicator on which any parallel communicator will occur.
//   context - The context used to compute the properties of the Jacobian J,
//             as well as any preconditioner.
//   cvode_type - A flag that determines the method to use for integration. 
//                Set to CV_ADAMS for the non-stiff Adams-Moulton method, 
//                and CV_BDF for the stiff backward difference formulae.
//   order - The desired order of accuracy of the integrator.
//   compute_Jv - A function that applies the Jacobian J to the vector x.
//   precond_type - A flag that determines the type of preconditioning to 
//                  be used. Possible values are PREC_NONE for no preconditioning,
//                  PREC_LEFT for left preconditioning, PREC_RIGHT for 
//                  right preconditioning, or PREC_BOTH for both left and right.
//   precond_setup - A function that sets up the preconditioner. NULL if precond_type is PREC_NONE.
//   precond_solve - A function that solves the preconditioner system. NULL if precond_type is PREC_NONE.
//   gram_schmidt - A flag that determines whether to use the classical 
//                  (CLASSICAL_GS) or the modified (MODIFIED_GS) Gram-Schmidt 
//                  orthogonalization process.
//   dtor - A destructor for freeing the context pointer (if any).
integrator_t* nonlinear_integrator_new(MPI_Comm comm,
                                       void* context, 
                                       int cvode_type,
                                       int order,
                                       integrator_compute_Jv_func compute_Jv,
                                       int precond_type,
                                       nonlinear_integrator_precond_setup_func precond_setup,
                                       nonlinear_integrator_precond_solve_func precond_solve,
                                       int gram_schmidt,
                                       integrator_dtor dtor);

#ifdef __cplusplus
}
#endif

#endif

