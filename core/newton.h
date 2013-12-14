// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

// This class solves (dense) systems of nonlinear equations using 
// Newton's method.
typedef struct newton_solver_t newton_solver_t;

// This function F(X) = 0 represents a nonlinear system of equations, and 
// uses the KINSOL interface.
typedef int (*newton_system_func_t)(N_Vector X, N_Vector F, void* context);

// Creates a new Newton solver for the given system of equations represented
// by the function F(X) = 0. The user-defined context is fed to the system 
// function when it is called. The last argument is an optional destructor 
// function for destroying the context when the Newton solver is destroyed.
newton_solver_t* newton_solver_new(int dimension,
                                   void* context,
                                   newton_system_func_t system_func,
                                   void (*context_dtor)(void*));

// Destroys the given Newton solver.
void newton_solver_free(newton_solver_t* solver);

// Returns the dimension of the system solved by the Newton solver.
int newton_solver_dimension(newton_solver_t* solver);

// Sets the tolerances for the function norm (norm_tolerance) and the Newton
// step (step_tolerance).
void newton_solver_set_tolerances(newton_solver_t* solver, double norm_tolerance, double step_tolerance);

// Sets the maximum number of Newton iterations for the solver.
void newton_solver_set_max_iterations(newton_solver_t* solver, int max_iterations);

// Given an initial guess, solve the system represented by the Newton solver.
// Returns true if the solve succeeded, false if not. num_iterations will 
// store the number of Newton iterations used to achieve the solution.
bool newton_solver_solve(newton_solver_t* solver, double* X, int* num_iterations);

// This variant of newton_solver_solve() lets one specify two vectors which can 
// scale the solution and the function to accelerate convergence:
// - x_scale: the diagonal components of a matrix Dx such that the components 
//            of Dx * x all have roughly the same magnitude as F(x) approaches 0.
// - F_scale: the diagonal components of a matrix Df such that the components 
//            of Df * F(x) all have roughly the same magnitude as F(x) approaches 0.
// If either of these arguments is NULL, the components of the corresponding 
// vector are assumed to be 1.
bool newton_solver_solve_scaled(newton_solver_t* solver, double* X, double* x_scale, double* F_scale, int* num_iterations);

// This defines a function that computes the value and derivative of a 
// (real-valued) function of a single variable given a context.
typedef void (*nonlinear_function_t)(void*, double, double*, double*);

// Solve the nonlinear equation F(x) = 0 using Newton-Raphson iteration on 
// the interval [min, max] using the initial estimate x, with the specified 
// maximum number of iterations and desired tolerance.
// Returns true if the solution was achieved, false otherwise.
bool newton_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int max_iters); 

// Solve the nonlinear equation F(x) = 0 using Brent's method on 
// the interval [x1, x2] with the specified maximum number of iterations 
// and desired tolerance. Returns the solution x. NOTE that the nonlinear
// function F need not have a derivative for Brent's method.
double brent_solve(nonlinear_function_t F, void* context, double x1, double x2, double tolerance, int max_iters, double* error); 

#endif

