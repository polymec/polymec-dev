// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KRYLOV_SOLVER_H
#define POLYMEC_KRYLOV_SOLVER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "core/sparse_local_matrix.h"

// Types of Krylov solver.
typedef enum
{
  KRYLOV_GMRES,
  KRYLOV_BICGSTAB,
  KRYLOV_TFQMR,
} krylov_t;

// Objects of this type solve linear systems A*x = b using preconditioned 
// Krylov subspace methods (built from matrix-vector products A*x, A*Ax, etc).
typedef struct krylov_solver_t krylov_solver_t;

// Creates a Krylov linear solver with a given maximum subspace dimension of 
// max_krylov_dim. matrix_vector_product is a function that computes the product 
// A*y for the linear operator A and the vector y, storing the result in the 
// array Ay. If the solver_type is KRYLOV_GMRES, the maximum number of restarts 
// is given by max_restarts--otherwise that parameter is ignored. The system has 
// num_local_values local equations, and num_remote_values remote ones.
krylov_solver_t* krylov_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   int (*matrix_vector_product)(void* context, real_t t, real_t* y, real_t* Ay),
                                   void (*dtor)(void* context),
                                   krylov_t solver_type,
                                   int max_krylov_dim,
                                   int max_restarts);

// Frees a solver.
void krylov_solver_free(krylov_solver_t* solver);

// Returns the context pointer for the solver.
void* krylov_solver_context(krylov_solver_t* solver);

// Returns the number of (local) equations in the linear system.
int krylov_solver_num_equations(krylov_solver_t* solver);

// Sets the tolerance for the residual norm (norm_tolerance) = |R|, where 
// R = A*x - b.
void krylov_solver_set_tolerance(krylov_solver_t* solver, real_t residual_tolerance);

// Sets up the Jacobi preconditioner for the Krylov solver, using the given 
// sparsity graph.
void krylov_solver_set_jacobi_preconditioner(krylov_solver_t* solver, 
                                             adj_graph_t* sparsity);

// Sets up a block Jacobi preconditioner for the Krylov solver with the given 
// block size and with the given sparsity graph.
void krylov_solver_set_block_jacobi_preconditioner(krylov_solver_t* solver, 
                                                   int block_size,
                                                   adj_graph_t* sparsity);

// Sets up an LU preconditioner for the Krylov solver with the given 
// sparsity graph.
void krylov_solver_set_lu_preconditioner(krylov_solver_t* solver, 
                                         adj_graph_t* sparsity);

// Sets up an incomplete LU preconditioner for the Krylov solver with the 
// given ILU parameters and sparsity graph.
void krylov_solver_set_ilu_preconditioner(krylov_solver_t* solver, 
                                          ilu_params_t* ilu_params,
                                          adj_graph_t* sparsity);

// Solves the linear system of equations A * x = b in place, storing the solution
// in the array b. Returns true if the solution was obtained, false if not. The 
// number of linear iterations will be stored in num_iterations upon success.
// If the solve is unsuccessful, the vector b will NOT contain the solution 
// unless the residual was reduced.
bool krylov_solver_solve(krylov_solver_t* solver, real_t t, real_t* b, 
                         real_t* residual_norm, int* num_iterations);

// Computes the residual of the linear system, R = A*x - b for the given value of x.
void krylov_solver_eval_residual(krylov_solver_t* solver, real_t t, real_t* x, real_t* b, real_t* R);

#endif

