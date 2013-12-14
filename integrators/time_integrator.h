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


#ifndef POLYMEC_TIME_INTEGRATOR_H
#define POLYMEC_TIME_INTEGRATOR_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "cvode/cvode.h"
#include "cvode/cvode_spils.h"

typedef enum
{
  GMRES,
  BICGSTAB
} time_integrator_solver_type_t;

typedef struct
{
  CVRhsFn rhs;
  void (*dtor)(void*);
  adj_graph_t* (*graph)(void*);
} time_integrator_vtable;

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. 
typedef struct time_integrator_t time_integrator_t;

// Creates a time integrator with the given name, context, and virtual table.
// Also given are some metadata, such as the order of the method, whether 
// it is explicit, implicit, or hybrid. The integrator will integrate 
// a system of N differential equations.
time_integrator_t* time_integrator_new(const char* name, 
                                       void* context,
                                       time_integrator_vtable vtable,
                                       int order,
                                       time_integrator_solver_type_t solver_type);

// Frees a time integrator.
void time_integrator_free(time_integrator_t* integrator);

// Returns the name of the integrator (internally stored).
char* time_integrator_name(time_integrator_t* integrator);

// Returns the context object for this integrator.
void* time_integrator_context(time_integrator_t* integrator);

// Returns the order of the integration method.
int time_integrator_order(time_integrator_t* integrator);

// Integrates the given solution X in place from time t1 to t2.
void time_integrator_step(time_integrator_t* integrator, double t1, double t2, double* X);

#endif

