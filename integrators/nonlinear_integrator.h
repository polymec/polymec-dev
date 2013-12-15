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

#ifndef POLYMEC_NONLINEAR_INTEGRATOR_H
#define POLYMEC_NONLINEAR_INTEGRATOR_H

#include "kinsol/kinsol.h"
#include "core/polymec.h"
#include "core/adj_graph.h"

// The different algorithms for matrix-free solution of nonlinear equations.
typedef enum
{
  GMRES,
  BICGSTAB,
  TFQMR
} nonlinear_integrator_type_t;

typedef struct
{
  KINSysFn eval;
  void (*dtor)(void*);
  adj_graph_t* (*graph)(void*);
} nonlinear_integrator_vtable;

// This class represents a collection of algorithms for integrating partial 
// differential equations that are discretized into a (sparse) system of 
// nonlinear equations. The integration is performed using Matrix-free 
// Newton-Krylov methods provided by KINSol.
typedef struct nonlinear_integrator_t nonlinear_integrator_t;

// Creates a solver with the given name, context, and virtual table for 
// solving a differential algebraic system with N equations.
nonlinear_integrator_t* nonlinear_integrator_new(const char* name,
                                                 void* context,
                                                 nonlinear_integrator_vtable vtable,
                                                 nonlinear_integrator_type_t type);

// Frees a solver.
void nonlinear_integrator_free(nonlinear_integrator_t* integrator);

// Returns an internal string storing the name of the solver.
char* nonlinear_integrator_name(nonlinear_integrator_t* integrator);

// Returns the context pointer for the solver.
void* nonlinear_integrator_context(nonlinear_integrator_t* integrator);

// Integrates the nonlinear system of equations F(X, t) = 0 in place, 
// using X as the initial guess.
void nonlinear_integrator_solve(nonlinear_integrator_t* integrator, double t, double* X);

#endif

