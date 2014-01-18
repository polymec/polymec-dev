// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_NEWTON_H
#define POLYMEC_NEWTON_H

#include "core/polymec.h"
#include "core/sundials_helpers.h"
#include "sundials/sundials_direct.h"

// Solves the (scalar) nonlinear equation F(x) = 0 using Brent's method on 
// the interval [x1, x2] with the specified maximum number of iterations 
// and desired tolerance. F is a function that takes a user-defined context and 
// a value x and returns F(x). This function returns the solution x which 
// satisfies fabs(F(x)) = tolerance.
real_t brent_solve(real_t (*F)(void*, real_t), void* context, real_t x1, real_t x2, real_t tolerance, int max_iters); 

// This class solves (dense) systems of nonlinear equations using 
// Newton's method.
typedef struct newton_solver_t newton_solver_t;

// This function F(X) = 0 represents a nonlinear system of equations, and 
// uses the KINSOL interface.
typedef int (*newton_system_func)(N_Vector X, N_Vector F, void* context);

// This function represents the Jacobian of a nonlinear system of equations, 
// and uses the KINSOL interface.
typedef int (*newton_jacobian_func)(long N, N_Vector X, N_Vector F, DlsMat J,
                                    void* context, N_Vector work1, N_Vector work2);

// Creates a new Newton solver for the given system of equations represented
// by the function F(X) = 0. The user-defined context is fed to the system 
// function when it is called. The last argument is an optional destructor 
// function for destroying the context when the Newton solver is destroyed.
newton_solver_t* newton_solver_new(int dimension,
                                   void* context,
                                   newton_system_func system_func,
                                   void (*context_dtor)(void*));

// This variant of newton_solver_new creates a Newton solver that uses the 
// additionally-supplied Jacobian function to accelerate the solution process.
newton_solver_t* newton_solver_new_with_jacobian(int dimension,
                                                 void* context,
                                                 newton_system_func system_func,
                                                 newton_jacobian_func jacobian_func,
                                                 void (*context_dtor)(void*));

// Destroys the given Newton solver.
void newton_solver_free(newton_solver_t* solver);

// Returns the dimension of the system solved by the Newton solver.
int newton_solver_dimension(newton_solver_t* solver);

// Sets the tolerances for the function norm (norm_tolerance) and the Newton
// step (step_tolerance).
void newton_solver_set_tolerances(newton_solver_t* solver, real_t norm_tolerance, real_t step_tolerance);

// Sets the maximum number of Newton iterations for the solver.
void newton_solver_set_max_iterations(newton_solver_t* solver, int max_iterations);

// Given an initial guess, solve the system represented by the Newton solver.
// Returns true if the solve succeeded, false if not. num_iterations will 
// store the number of Newton iterations used to achieve the solution.
bool newton_solver_solve(newton_solver_t* solver, real_t* X, int* num_iterations);

// This variant of newton_solver_solve() lets one specify two vectors which can 
// scale the solution and the function to accelerate convergence:
// - x_scale: the diagonal components of a matrix Dx such that the components 
//            of Dx * x all have roughly the same magnitude as F(x) approaches 0.
// - F_scale: the diagonal components of a matrix Df such that the components 
//            of Df * F(x) all have roughly the same magnitude as F(x) approaches 0.
// If either of these arguments is NULL, the components of the corresponding 
// vector are assumed to be 1.
bool newton_solver_solve_scaled(newton_solver_t* solver, real_t* X, real_t* x_scale, real_t* F_scale, int* num_iterations);

#endif

