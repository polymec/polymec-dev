// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_KRYLOV_SPARSE_LIN_SOLVERS_H
#define POLYMEC_KRYLOV_SPARSE_LIN_SOLVERS_H

#include "core/sparse_lin_solver.h"
#include "core/sundials_helpers.h"

// A prototype for a function that computes a matrix-vector product for a
// linear system. This should return 0 on success and non-zero on failure.
typedef int (*mf_krylov_sparse_lin_solver_compute_Ax_func)(void *A_data, N_Vector x, N_Vector Ax);

// A prototype for a function that solves the preconditioner system
// M * z = r. Here, precond_type = PRECOND_LEFT for left-preconditioned 
// systems and PRECOND_RIGHT for right-preconditioned systems. This should 
// return 0 on success and non-zero on failure.
typedef int (*mf_krylov_sparse_lin_solver_precond_func)(void *P_data, N_Vector r, N_Vector z, int precond_type);

// A prototype for a function that computes scaling factors s1 and s2 in 
// the preconditioned linear system.
typedef void (*mf_krylov_sparse_lin_solver_scaling_func)(void *context, N_Vector s1, N_Vector s2);

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
//   compute_scale_factors - A function that computes scaling factors s1, s2.
//                           These can be NULL if no scaling is required.
//   delta - The tolerance on the L2 norm of the scaled, preconditioned 
//           residual.
//   max_restarts - The maximum number of times the algorithm is allowed 
//                  to restart.
//   dtor - A destructor for freeing the context pointer (if any).
sparse_lin_solver_t* mf_gmres_sparse_lin_solver_new(MPI_Comm comm,
                                                    void* context, 
                                                    mf_krylov_sparse_lin_solver_compute_Ax_func Ax,
                                                    int max_kdim,
                                                    int gram_schmidt,
                                                    int precond_type,
                                                    mf_krylov_sparse_lin_solver_precond_func precond,
                                                    mf_krylov_sparse_lin_solver_scaling_func compute_scale_factors,
                                                    double delta,
                                                    int max_restarts,
                                                    sparse_lin_solver_dtor dtor);

// Creates a matrix-free sparse linear solver that uses the stabilized Bi-conjugate
// Gradient (Bi_CGSTAB) method. Arguments:
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
//   compute_scale_factors - A function that computes scaling factors s1, s2.
//                           These can be NULL if no scaling is required.
//   delta - The tolerance on the L2 norm of the scaled, preconditioned 
//           residual.
//   dtor - A destructor for freeing the context pointer (if any).
sparse_lin_solver_t* mf_bicgstab_sparse_lin_solver_new(MPI_Comm comm,
                                                       void* context, 
                                                       mf_krylov_sparse_lin_solver_compute_Ax_func Ax,
                                                       int max_kdim,
                                                       int precond_type,
                                                       mf_krylov_sparse_lin_solver_precond_func precond,
                                                       mf_krylov_sparse_lin_solver_scaling_func compute_scale_factors,
                                                       double delta,
                                                       sparse_lin_solver_dtor dtor);

#endif

