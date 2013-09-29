// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_INTEGRATOR_H
#define POLYMEC_INTEGRATOR_H

#include "core/polymec.h"

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. 
typedef struct integrator_t integrator_t;

// A prototype for a function that computes the right hand side RHS of the ith
// differential equation within a nonlinear system of N equations at time t.
typedef void (*integrator_compute_rhs_func)(void* context, int i, double t, double* X, int N, double* RHS);

// A function for initializing the solution X to a nonlinear system of 
// N differential equations at a given time t.
typedef void (*integrator_init_func)(void* context, double t, int N, double* X0);

// A function for integrating a given solution vector for a system of N 
// nonlinear equations from time t1 to t2.
typedef void (*integrator_step_func)(void* context, double t1, double t2, int N, double* X);

// A destructor function for the context object (if any).
typedef void (*integrator_dtor)(void*);

// This virtual table must be implemented by any integrator.
typedef struct 
{
  integrator_init_func        init;
  integrator_compute_rhs_func compute_rhs;
  integrator_step_func        step;
  integrator_dtor             dtor;
} integrator_vtable;

// Types of integrators (for algorithmic classification).
typedef enum
{
  INTEGRATOR_EXPLICIT,
  INTEGRATOR_IMPLICIT,
  INTEGRATOR_SEMI_IMPLICIT
} integrator_type_t;

// Creates an integrator with the given name, context, and virtual table.
// Also given are some metadata, such as the order of the method, whether 
// it is explicit, implicit, or hybrid. The integrator will integrate 
// a system of N differential equations.
integrator_t* integrator_new(const char* name, 
                             void* context,
                             integrator_vtable vtable,
                             int order,
                             integrator_type_t type, 
                             int N);

// Frees an integrator.
void integrator_free(integrator_t* integrator);

// Returns the name of the integrator (internally stored).
char* integrator_name(integrator_t* integrator);

// Returns the context object for this integrator.
void* integrator_context(integrator_t* integrator);

// Returns the order of the integration method.
int integrator_order(integrator_t* integrator);

// Returns the type of the integrator (explicit, implicit, IMEX).
integrator_type_t integrator_type(integrator_t* integrator);

// Returns the number of differential equations integrated by the integrator.
int integrator_N(integrator_t* integrator);

// Initializes the integrator and the solution vector X0 at time t.
void integrator_init(integrator_t* integrator, double t, double* X0);

// Integrates the given solution X in place from time t1 to t2.
void integrator_step(integrator_t* integrator, double t1, double t2, double* X);

#endif

