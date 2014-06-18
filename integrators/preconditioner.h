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

#ifndef POLYMEC_PRECONDITIONER_H
#define POLYMEC_PRECONDITIONER_H

#include "core/polymec.h"
#include "core/adj_graph.h"

// This class constructs and solves a preconditioned linear system to aid 
// the Jacobian-free Newton-Krylov solution of a nonlinear system.
typedef struct preconditioner_t preconditioner_t;

// This function sets up the preconditioner to obtain a preconditioned solution.
typedef void (*preconditioner_setup_func)(void* context, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* xdot);

// This function solves the preconditioner system A*X = B, where A is the 
// preconditioner matrix and B is the given right-hand side. The right-hand 
// side is taken as input and is filled with the solution on output.
typedef bool (*preconditioner_solve_func)(void* context, real_t* B);

// This function writes a text representation of the preconditioner to the given file.
typedef void (*preconditioner_fprintf_func)(void* context, FILE* stream);

// This function destroys the data context for the preconditioner.
typedef void (*preconditioner_dtor)(void* context);

// This virtual table determines the behavior of the preconditioner.
typedef struct
{
  preconditioner_setup_func   setup;
  preconditioner_solve_func   solve;
  preconditioner_fprintf_func fprintf;
  preconditioner_dtor         dtor;
} preconditioner_vtable;

// Constructs a preconditioner with the given name, data context, 
// virtual table. 
preconditioner_t* preconditioner_new(const char* name,
                                     void* context,
                                     preconditioner_vtable vtable);

// Destroys the given preconditioner.
void preconditioner_free(preconditioner_t* precond);

// Returns an internal pointer to the name of the preconditioner.
char* preconditioner_name(preconditioner_t* precond);

// Returns the data context for the preconditioner.
void* preconditioner_context(preconditioner_t* precond);

// Sets up the preconditioner to obtain a representation of the preconditioner
// operator. This operator represents the expression 
//
// alpha I + beta * dF/dx + gamma * dF/dxdot
//
// where I is the identity matrix, and dF/dx and dF/dxdot are the derivatives 
// of a system function F with respect to the solution and its time derivative.
// Note that xdot can be NULL if gamma == 0.0.
void preconditioner_setup(preconditioner_t* precond, real_t alpha, real_t beta, real_t gamma,
                          real_t, real_t* x, real_t* xdot);

// Solves the preconditioned linear system. On input, rhs contains the 
// right-hand side of the system, and on output, it contains the solution 
// to the preconditioned system. Returns true if the system was successfully 
// solved, false if not.
bool preconditioner_solve(preconditioner_t* precond,
                          real_t* rhs);

// Writes a text representation of the given preconditioner to the given file.
void preconditioner_fprintf(preconditioner_t* precond, FILE* stream);

#endif

