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

#ifndef POLYMEC_nonlinear_SOLVER_H
#define POLYMEC_nonlinear_SOLVER_H

#include "core/polymec.h"

// A function for computing the value of the ith equation within an 
// N-dimensional differential algebraic system F(X) = 0 given a solution X at time t.
typedef double (*nonlinear_solver_eval_func)(void* context, 
                                             double t, 
                                             int i, 
                                             int N,
                                             double* X);

// A function for solving a system of differential algebraic equations at 
// a time t (in place).
typedef void (*nonlinear_solver_solve_func)(void* context,
                                            double t, 
                                            int N,
                                            double* X);

// A destructor for nonlinear solvers.
typedef void (*nonlinear_solver_dtor)(void* context);

// This virtual table must be implemented by any nonlinear_solver_t.
typedef struct 
{
  nonlinear_solver_eval_func  eval;
  nonlinear_solver_solve_func solve;
  nonlinear_solver_dtor       dtor;
} nonlinear_solver_vtable;

// This class represents a way of integrating a system of 
// nonlinear equations.
typedef struct nonlinear_solver_t nonlinear_solver_t;

// Creates a solver with the given name, context, and virtual table for 
// solving a differential algebraic system with N equations.
nonlinear_solver_t* nonlinear_solver_new(const char* name,
                             void* context,
                             nonlinear_solver_vtable vtable,
                             int N);

// Frees a solver.
void nonlinear_solver_free(nonlinear_solver_t* solver);

// Returns an internal string storing the name of the solver.
char* nonlinear_solver_name(nonlinear_solver_t* solver);

// Returns the context pointer for the solver.
void* nonlinear_solver_context(nonlinear_solver_t* solver);

// Solves the system of equations F(X, t) = 0 in place, using X as the initial guess.
void nonlinear_solver_solve(nonlinear_solver_t* solver, double t, double* X);

#endif

