#ifndef POLYMEC_BDF_INTEGRATOR_H
#define POLYMEC_BDF_INTEGRATOR_H

#include "integrators/nonlinear_integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a integrator that integrates a nonlinear system of equations using 
// the stiff backward difference formulae implemented in CVODE.
//   comm - The MPI communicator on which any parallel communicator will occur.
//   context - The context used to compute the properties of the Jacobian J,
//             as well as any preconditioner.
//   order - The desired order of accuracy of the integrator.
//   compute_F - A function that computes the nonlinear function F.
//   compute_Jv - A function that applies the Jacobian J to the vector x.
//   solver_type - A flag that indicates what kind of linear solver to use 
//                 (INTEGRATOR_SOLVER_GMRES or INTEGRATOR_SOLVER_BICGSTAB).
//   max_kdim - The maximum dimension of the Krylov space used by the solver.
//   gram_schmidt - A flag that determines whether to use the classical 
//                  (CLASSICAL_GS) or the modified (MODIFIED_GS) Gram-Schmidt 
//                  orthogonalization process.
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
integrator_t* bdf_integrator_new(MPI_Comm comm,
                                 void* context, 
                                 int order,
                                 integrator_compute_F_func compute_F,
                                 integrator_compute_Jv_func compute_Jv,
                                 nonlinear_integrator_solver_type_t solver_type,
                                 int max_kdim,
                                 int gram_schmidt,
                                 int precond_type,
                                 nonlinear_integrator_precond_setup_func precond_setup,
                                 nonlinear_integrator_precond_solve_func precond_solve,
                                 integrator_dtor dtor);

#ifdef __cplusplus
}
#endif

#endif

