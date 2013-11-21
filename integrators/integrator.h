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
#include "core/adj_graph.h"
#include "cvode/cvode.h"
#include "cvode/cvode_spils.h"

typedef enum
{
  GMRES,
  BICGSTAB
} integrator_solver_type_t;

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. 
typedef struct integrator_t integrator_t;

// Creates an integrator with the given name, context, and virtual table.
// Also given are some metadata, such as the order of the method, whether 
// it is explicit, implicit, or hybrid. The integrator will integrate 
// a system of N differential equations.
integrator_t* integrator_new(const char* name, 
                             void* context,
                             CVRhsFn rhs,
                             void (*dtor)(void*),
                             adj_graph_t* graph,
                             int order,
                             integrator_solver_type_t solver_type);

// Frees an integrator.
void integrator_free(integrator_t* integrator);

// Returns the name of the integrator (internally stored).
char* integrator_name(integrator_t* integrator);

// Returns the context object for this integrator.
void* integrator_context(integrator_t* integrator);

// Returns the order of the integration method.
int integrator_order(integrator_t* integrator);

// Integrates the given solution X in place from time t1 to t2.
void integrator_step(integrator_t* integrator, double t1, double t2, double* X);

#endif

