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

// This defines a vector-valued function that computes its value 
// given a context. 
typedef void (*nonlinear_vector_function_t)(void*, double*, double*);

// This function prototype describes functions that compute Jacobians.
typedef void (*jacobian_function_t)(void*, nonlinear_vector_function_t, int, double*, double*);

// This struct represents a nonlinear system of equations F(x) = 0.
typedef struct
{
  int dim;                                // Dimension of the system.
  nonlinear_vector_function_t compute_F;  // Vector-valued function F.
  jacobian_function_t compute_J;          // Jacobian calculation (optional).
  void* context;                          // Context pointer.
} nonlinear_system_t;

// Solve the nonlinear system F(x) = 0 using a globally-convergent version 
// of Newton's method with the initial estimate vector x, with the specified 
// maximum number of iterations and desired tolerance. Returns true if the 
// solution converged, false otherwise. If check_solution is set to true, it is possible that spurious convergence occurred.
bool newton_solve_system(nonlinear_system_t* system, double* x, double tolerance, int max_iters, int* num_iters);

// Solve the nonlinear system F(x) = 0 using Broyden's method with the 
// initial estimate vector x, with the specified maximum number of iterations 
// and desired tolerance. Returns true if the solution converged, false otherwise. 
// If check_solution is set to true, it is possible that spurious convergence occurred.
bool broyden_solve_system(nonlinear_system_t* system, double* x, double tolerance, int max_iters, int* num_iters);

#endif

