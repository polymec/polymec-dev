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

#ifndef POLYMEC_NEWTON_PRECONDITIONER_H
#define POLYMEC_NEWTON_PRECONDITIONER_H

#include "core/polymec.h"
#include "core/preconditioner.h"

// The machinery in this file is used to construct preconditioners for 
// large distributed Newton methods of the following types:
// (1) Finding the solution x of a nonlinear system F(t, x) = 0.
// (2) Integrating a nonlinear system of ODEs dx/dt = F(t, x).
// (3) Solving a system of nonlinear differential algebraic equations (DAE):
//     F(t, x, xdot) = 0.
// Newton preconditioners are subclasses of the basic preconditioner class.

// Every Newton preconditioner has to implement the functions in the 
// following virtual table.
typedef struct
{
  // Method for computing the preconditioner operator 
  // P = alpha * I + beta * dF/dx + gamma * dF/dxdot.
  // where F(t, x, xdot) is a nonlinear function whose solution x is sought.
  void (*compute_P)(void* context, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* xdot);
  // Method for solving the preconditioner system P * z = r. On input, z is 
  // the right-hand side, and on output it is the solution.
  bool (*solve)(void* context, real_t* z);
  // fprintf function.
  void (*fprintf)(void* context, FILE* stream);
  // Context destructor.
  void (*dtor)(void* context);
} newton_preconditioner_vtable;

// Constructs a preconditioner that represents the linear combination
// alpha * I + beta * dF/dx + gamma * dF/dx, where F(t, x[, xdot]) = 0 is a 
// function representing a nonlinear ODE [or differential/algebraic) system, 
// I is the identity matrix, and dF/dx [and dF/dxdot] are the derivatives of 
// F with respect to x [and xdot].
preconditioner_t* newton_preconditioner_new(const char* name,
                                            void* context,
                                            newton_preconditioner_vtable vtable);

// Sets up the preconditioner at the point (t, x, xdot) in solution space, computing
// alpha * I + beta * dF/dx + gamma * dF/dx for the given alpha, beta, and gamma.
// Note that xdot can be NULL if F is only a function of x.
void newton_preconditioner_setup(preconditioner_t* precond, 
                                 real_t alpha, real_t beta, real_t gamma,
                                 real_t t, real_t* x, real_t* xdot);

#endif

