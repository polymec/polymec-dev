// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KRYLOV_SOLVER_H
#define POLYMEC_KRYLOV_SOLVER_H

#include "core/polymec.h"
#include "core/preconditioner.h"

// A Krylov solver is an iterative linear solver that uses Krylov subspace 
// methods to solve a linear system A * X = B.
typedef struct krylov_solver_t krylov_solver_t;

// This virtual table must be defined by any implementation of a Krylov solver.
typedef struct
{
  int (*ax)(void* context, real_t* x, real_t* Ax, int N); // Implements A*X.
  void (*b)(void* context, real_t* b, int N); // Computes B, the right hand side.
  void (*dtor)(void* context); // context destructor
} krylov_solver_vtable;

// Creates a Krylov solver that solves A * X = B using a Generalized Minimum
// Residual (GMRES) method with the given maximum subspace dimension, the given 
// maximum number of restarts, and the given problem size N. The given context 
// and virtual table define the behavior of the solver. 
krylov_solver_t* gmres_krylov_solver_new(MPI_Comm comm,
                                         void* context,
                                         krylov_solver_vtable vtable,
                                         int N,
                                         int max_krylov_dim,
                                         int max_restarts);
  
// Creates a Krylov solver that solves A * X = B using a Stabilized Bi-Conjugate
// Gradient (BiCGSTAB) method with the given maximum subspace dimension and 
// the given problem size N. The given context and virtual table define the 
// behavior of the solver.
krylov_solver_t* bicgstab_krylov_solver_new(MPI_Comm comm,
                                            void* context,
                                            krylov_solver_vtable vtable,
                                            int N,
                                            int max_krylov_dim);
  
// Creates a Krylov solver that solves A * X = B using a Transpose­Free 
// Quasi-Minimum Residual (TFQMR) method with the given maximum subspace 
// dimension and the given problem size N. The given context and virtual table 
// define the behavior of the solver.
krylov_solver_t* tfqmr_krylov_solver_new(MPI_Comm comm,
                                         void* context,
                                         krylov_solver_vtable vtable,
                                         int N,
                                         int max_krylov_dim);
 
// Frees the resources associated with the given Krylov solver.
void krylov_solver_free(krylov_solver_t* solver);

// Returns the context pointer associated with this solver.
void* krylov_solver_context(krylov_solver_t* solver);

// Sets the tolerance on the L2 norm of the residual for the Krylov solver.
// By default, the tolerance is 1e-8.
void krylov_solver_set_tolerance(krylov_solver_t* solver, real_t delta);

// Sets the (linear) preconditioner to use for this Krylov solver.
void krylov_solver_set_preconditioner(krylov_solver_t* solver, preconditioner_t* precond);

// Returns an internal pointer to the diagonal scaling matrix S1. By default,
// S2 is the identity matrix.
real_t* krylov_solver_s1(krylov_solver_t* solver);

// Returns an internal pointer to the diagonal scaling matrix S2. By default,
// S2 is the identity matrix.
real_t* krylov_solver_s2(krylov_solver_t* solver);

// Returns a newly-created array that is large enough to store the solution 
// vector for the linear system represented by the Krylov solver. Contains 
// zeros. Must be freed using polymec_free().
real_t* krylov_solver_vector(krylov_solver_t* solver);

// Solves the linear system A * X = B defined by the Krylov solver, storing:
// - the solution in X,
// - the L2 norm of the residual in res_norm,
// - the number of iterations in num_iters.
// - the number of preconditioning calls in num_precond.
// Returns true if the linear system was solved (however approximately), 
// false if not.
bool krylov_solver_solve(krylov_solver_t* solver, real_t* X, real_t* res_norm, 
                         int* num_iters, int* num_precond);

#endif
