// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NEWTON_H
#define POLYMEC_NEWTON_H

#include "core/polymec.h"

// Solves the (scalar) nonlinear equation F(x) = 0 using Brent's method on 
// the interval [x1, x2] with the specified maximum number of iterations 
// and desired tolerance. F is a function that takes a user-defined context and 
// a value x and returns F(x). This function returns the solution x which 
// satisfies ABS(F(x)) = tolerance.
real_t brent_solve(real_t (*F)(void*, real_t), void* context, real_t x1, real_t x2, real_t tolerance, int max_iters); 

// This class solves (dense) systems of nonlinear equations using 
// Newton's method.
typedef struct dense_newton_solver_t dense_newton_solver_t;

// This function F(X) = 0 represents a nonlinear system of equations.
// It should return 0 on success, 1 on a recoverable error, and -1 on a 
// fatal error.
typedef int (*dense_newton_system_func)(void* context, real_t* X, real_t* F);

// This function represents the Jacobian of a nonlinear system of equations.
// N is the dimension of the Jacobian.
// X is the solution about which the system is linearized.
// F contains the values of the system function at X.
// work1 and work2 are work vectors (sized to match X and F).
// J is a matrix stored in column-major order within an N*N array.
// This function should return 0 on success, nonzero otherwise.
typedef int (*dense_newton_jacobian_func)(void* context, int N, real_t* X, real_t* F, 
                                          real_t* work1, real_t* work2, real_t* J);

// Creates a new dense Newton solver for the given system of equations represented
// by the function F(X) = 0. The user-defined context is fed to the system 
// function when it is called. The last argument is an optional destructor 
// function for destroying the context when the Newton solver is destroyed.
dense_newton_solver_t* dense_newton_solver_new(int dimension,
                                               void* context,
                                               dense_newton_system_func system_func,
                                               void (*context_dtor)(void*));

// This variant of dense_newton_solver_new creates a Newton solver that uses the 
// additionally-supplied Jacobian function to accelerate the solution process.
dense_newton_solver_t* dense_newton_solver_new_with_jacobian(int dimension,
                                                             void* context,
                                                             dense_newton_system_func system_func,
                                                             dense_newton_jacobian_func jacobian_func,
                                                             void (*context_dtor)(void*));

// Destroys the given dense Newton solver.
void dense_newton_solver_free(dense_newton_solver_t* solver);

// Returns the dimension of the system solved by the dense Newton solver.
int dense_newton_solver_dimension(dense_newton_solver_t* solver);

// Sets the tolerances for the function norm (norm_tolerance) and the Newton
// step (step_tolerance).
void dense_newton_solver_set_tolerances(dense_newton_solver_t* solver, real_t norm_tolerance, real_t step_tolerance);

// Sets the maximum number of Newton iterations for the solver.
void dense_newton_solver_set_max_iterations(dense_newton_solver_t* solver, int max_iterations);

// Given an initial guess, solve the system represented by the Newton solver.
// Returns true if the solve succeeded, false if not. num_iterations will 
// store the number of Newton iterations used to achieve the solution.
bool dense_newton_solver_solve(dense_newton_solver_t* solver, real_t* X, int* num_iterations);

// This variant of newton_solver_solve() lets one specify two vectors which can 
// scale the solution and the function to accelerate convergence:
// - x_scale: the diagonal components of a matrix Dx such that the components 
//            of Dx * x all have roughly the same magnitude as F(x) approaches 0.
// - F_scale: the diagonal components of a matrix Df such that the components 
//            of Df * F(x) all have roughly the same magnitude as F(x) approaches 0.
// If either of these arguments is NULL, the components of the corresponding 
// vector are assumed to be 1.
bool dense_newton_solver_solve_scaled(dense_newton_solver_t* solver, real_t* X, real_t* x_scale, real_t* F_scale, int* num_iterations);

#endif

