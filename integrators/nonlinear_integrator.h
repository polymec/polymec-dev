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

#ifndef POLYMEC_NONLINEAR_INTEGRATOR_H
#define POLYMEC_NONLINEAR_INTEGRATOR_H

#include "core/polymec.h"
#include "core/adj_graph.h"

// The different algorithms for matrix-free solution of nonlinear equations.
typedef enum
{
  GMRES,
  BICGSTAB,
  TFQMR
} nonlinear_integrator_type_t;

// This virtual table determines the behavior of the nonlinear integrator.
typedef struct
{
  // This function evaluates the residual function for the nonlinear system
  // of equations using the solution vector x and placing the result in F.
  int (*eval)(void* context, real_t* x, real_t* F);

  // This (optional) function allows the state (context) to set the time at 
  // which the equations are to be integrated.
  void (*set_time)(void* context, real_t t);

  // This (optional) function sets the "x-scaling vector," which contains the diagonal 
  // components of a matrix Dx such that the components of Dx * x all have 
  // roughly the same magnitude as F(x) approaches 0.
  // - F_scale: the diagonal components of a matrix Df such that the components 
  //            of Df * F(x) all have roughly the same magnitude as F(x) approaches 0.
  void (*set_x_scale)(void* context, real_t* x_scale);

  // This (optional) function sets the "F-scaling vector," which contains the diagonal 
  // components of a matrix Df such that the components of Df * F(x) all have 
  // roughly the same magnitude as F(x) approaches 0.
  void (*set_F_scale)(void* context, real_t* F_scale);

  // This (optional) function sets the contraints vector, which places algebraic 
  // constraints on the components of the solution vector x. If constraints[i] is:
  // 0.0  - no constraint is placed on x[i].
  // 1.0  - x[i] must be non-negative.
  // -1.0 - x[i] must be non-positive.
  // 2.0  - x[i] must be positive.
  void (*set_constraints)(void* context, real_t* constraints);

  // This (optional) function destroys the state (context) when the nonlinear integrator 
  // is destroyed.
  void (*dtor)(void* context);

  // This function returns the adjacency graph reflecting the sparsity of the 
  // nonlinear system. It is *not* a block graph, so any nonzero blocks should 
  // be reflected as groups of vertices in the graph.
  adj_graph_t* (*graph)(void* context);

} nonlinear_integrator_vtable;

// This class represents a collection of algorithms for integrating partial 
// differential equations that are discretized into a (sparse) system of 
// nonlinear equations. The integration is performed using Matrix-free 
// Newton-Krylov methods provided by KINSol.
typedef struct nonlinear_integrator_t nonlinear_integrator_t;

// Creates an integrator with the given name, context, and virtual table for 
// integrated a discretized system of partial differential equations
// represented by a sparse nonlinear system of equations. The virtual table 
// defines accessor methods for the residual function and the adjacency graph, 
// and type indicates which type of Krylov method is used to solve the underlying 
// linear systems (using a maximum subspace dimension of max_krylov_dim).
nonlinear_integrator_t* nonlinear_integrator_new(const char* name,
                                                 void* context,
                                                 MPI_Comm comm,
                                                 nonlinear_integrator_vtable vtable,
                                                 nonlinear_integrator_type_t type,
                                                 int max_krylov_dim);

// Frees a solver.
void nonlinear_integrator_free(nonlinear_integrator_t* integrator);

// Returns an internal string storing the name of the solver.
char* nonlinear_integrator_name(nonlinear_integrator_t* integrator);

// Returns the context pointer for the solver.
void* nonlinear_integrator_context(nonlinear_integrator_t* integrator);

// Sets the tolerances for the function norm (norm_tolerance) and the Newton
// step (step_tolerance) for the nonlinear integrator.
void newton_solver_set_tolerances(nonlinear_integrator_t* integrator, real_t norm_tolerance, real_t step_tolerance);

// Sets the maximum number of Newton iterations for the integrator.
void newton_solver_set_max_iterations(nonlinear_integrator_t* integrator, int max_iterations);

// Integrates the nonlinear system of equations F(X, t) = 0 in place, 
// using X as the initial guess. Returns true if the solution was obtained, 
// false if not. The number of nonlinear iterations will be stored in 
// num_iterations upon success.
bool nonlinear_integrator_solve(nonlinear_integrator_t* integrator, real_t t, real_t* X, int* num_iterations);

#endif

