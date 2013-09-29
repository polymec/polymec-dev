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

#ifndef POLYMEC_DAE_SOLVER_H
#define POLYMEC_DAE_SOLVER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "core/table.h"

// A function for computing the ith equation within an N-dimensional 
// differential algebraic system F(X) = 0 given a solution X at time t.
typedef void (*dae_solver_eval_dae_func)(void* context, 
                                         double t, 
                                         int N,
                                         double* X, 
                                         double* F);

// A function for solving a system of differential algebraic equations at 
// a time t (in place).
typedef void (*dae_solver_solve_func)(void* context,
                                      double t, 
                                      int N,
                                      double* X);

// A destructor for DAE solvers.
typedef void (*dae_solver_dtor)(void* context);

// This virtual table must be implemented by any dae_solver_t.
typedef struct 
{
  dae_solver_eval_dae_func eval_dae;
  dae_solver_solve_func    solve;
  dae_solver_dtor          dtor;
} dae_solver_vtable;

// This class represents a way of integrating a system of 
// differential algebraic equations (DAE).
typedef struct dae_solver_t dae_solver_t;

// Creates a solver with the given name, context, and virtual table for 
// solving a differential algebraic system with N equations.
dae_solver_t* dae_solver_new(const char* name,
                             void* context,
                             dae_solver_vtable vtable,
                             int N);

// Frees a DAE solver.
void dae_solver_free(dae_solver_t* solver);

// Returns an internal string storing the name of the solver.
char* dae_solver_name(dae_solver_t* solver);

// Returns the context pointer for the solver.
void* dae_solver_context(dae_solver_t* solver);

// Solves the system of equations F(X, t) = 0 in place, using X as the initial guess.
void dae_solver_solve(dae_solver_t* solver, double t, double* X);

#endif

