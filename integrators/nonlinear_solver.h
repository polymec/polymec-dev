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

#ifndef POLYMEC_NONLINEAR_SOLVER_H
#define POLYMEC_NONLINEAR_SOLVER_H

#include "kinsol/kinsol.h"
#include "core/polymec.h"
#include "core/adj_graph.h"

// The different algorithms for matrix-free solution of nonlinear equations.
typedef enum
{
  GMRES,
  BICGSTAB,
  TFQMR
} nonlinear_solver_type_t;

typedef struct
{
  KINSysFn eval;
  void (*dtor)(void*);
  adj_graph_t* (*graph)(void*);
} nonlinear_solver_vtable;

// This class represents a way of integrating a system of 
// nonlinear equations.
typedef struct nonlinear_solver_t nonlinear_solver_t;

// Creates a solver with the given name, context, and virtual table for 
// solving a differential algebraic system with N equations.
nonlinear_solver_t* nonlinear_solver_new(const char* name,
                                         void* context,
                                         nonlinear_solver_vtable vtable,
                                         nonlinear_solver_type_t type);

// Frees a solver.
void nonlinear_solver_free(nonlinear_solver_t* solver);

// Returns an internal string storing the name of the solver.
char* nonlinear_solver_name(nonlinear_solver_t* solver);

// Returns the context pointer for the solver.
void* nonlinear_solver_context(nonlinear_solver_t* solver);

// Solves the system of equations F(X, t) = 0 in place, using X as the initial guess.
void nonlinear_solver_solve(nonlinear_solver_t* solver, double t, double* X);

#endif

