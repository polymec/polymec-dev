#ifndef POLYMEC_KRYLOV_SPARSE_LIN_SOLVERS_H
#define POLYMEC_KRYLOV_SPARSE_LIN_SOLVERS_H

#include "core/sparse_lin_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

// A prototype for a function that computes a matrix-vector product for a
// linear system.
typedef int (*krylov_sparse_lin_solver_compute_Ax_func)(void *A_data, N_Vector x, N_Vector Ax);

// A prototype for a function that solves the preconditioner system
// M * z = r. Here, precond_type = PRECOND_LEFT for left-preconditioned 
// systems and PRECOND_RIGHT for right-preconditioned systems.
typedef int (*krylov_sparse_lin_solver_precond_func)(void *P_data, N_Vector r, N_Vector z, int precond_type);

// Creates a sparse linear solver that uses the Generalized Minimum 
// Residual (GMRES) method. Arguments:
//   comm - The MPI communicator on which any parallel communicator will occur.
//   context - The context used to compute the properties of the linear 
//             operator A and the right-hand-side vector b, as well as 
//             any preconditioner.
//   compute_Ax - A function that applies A to the vector x.
//   max_kdim - The maximum Krylov dimension that the integrator can use.
//   gram_schmidt - A flag that determines whether to use the classical 
//                  (CLASSICAL_GS) or the modified (MODIFIED_GS) Gram-Schmidt 
//                  orthogonalization process.
//   precond_type - A flag that determines the type of preconditioning to 
//                  be used. Possible values are PREC_NONE for no preconditioning,
//                  PREC_LEFT for left preconditioning, PREC_RIGHT for 
//                  right preconditioning, or PREC_BOTH for both left and right.
//   precond - A function that solves the preconditioner system Pz = r. 
//             This can be NULL if precond_type is PREC_NONE.
//   delta - The tolerance on the L2 norm of the scaled, preconditioned 
//           residual.
//   max_restarts - The maximum number of times the algorithm is allowed 
//                  to restart.
//   dtor - A destructor for freeing the context pointer (if any).
sparse_lin_solver_t* gmres_sparse_lin_solver_new(MPI_Comm comm,
                                                 void* context, 
                                                 krylov_sparse_lin_solver_compute_Ax_func Ax,
                                                 int max_kdim,
                                                 int gram_schmidt,
                                                 int precond_type,
                                                 krylov_sparse_lin_solver_precond_func precond,
                                                 double delta,
                                                 int max_restarts,
                                                 sparse_lin_solver_dtor dtor);

// Creates an integrator that uses the 1st-order implicit backward Euler 
// method to integrate a linear system of ordinary differential equations.
// This version solves the linear system using the stabilized Bi-conjugate
// Gradient (Bi-CGSTAB) method. Arguments:
//   comm - The MPI communicator on which any parallel communicator will occur.
//   context - The context used to compute the properties of the linear 
//             operator A and the right-hand-side vector b, as well as 
//             any preconditioner.
//   compute_Ax - A function that applies A to the vector x.
//   max_kdim - The maximum Krylov dimension that the integrator can use.
//   precond_type - A flag that determines the type of preconditioning to 
//                  be used. Possible values are PREC_NONE for no preconditioning,
//                  PREC_LEFT for left preconditioning, PREC_RIGHT for 
//                  right preconditioning, or PREC_BOTH for both left and right.
//   precond - A function that solves the preconditioner system Pz = r. 
//             This can be NULL if precond_type is PREC_NONE.
//   delta - The tolerance on the L2 norm of the scaled, preconditioned 
//           residual.
//   dtor - A destructor for freeing the context pointer (if any).
sparse_lin_solver_t* bicgstab_sparse_lin_solver_new(MPI_Comm comm,
                                                    void* context, 
                                                    krylov_sparse_lin_solver_compute_Ax_func Ax,
                                                    int max_kdim,
                                                    int precond_type,
                                                    krylov_sparse_lin_solver_precond_func precond,
                                                    double delta,
                                                    sparse_lin_solver_dtor dtor);

#ifdef __cplusplus
}
#endif

#endif

