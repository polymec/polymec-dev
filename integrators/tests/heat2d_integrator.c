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

//------------------------------------------------------------------------
//                        2D heat equation problem
//------------------------------------------------------------------------
// This is adapted from the IDA heat 2D example problem programmed by
// Allan Taylor, Alan Hindmarsh and Radu Serban @ LLNL
// This example solves a discretized 2D heat equation problem.
// This version uses the Krylov solver IDASpgmr.
//
// The DAE system solved is a spatial discretization of the PDE
//          du/dt = d^2u/dx^2 + d^2u/dy^2
// on the unit square. The boundary condition is u = 0 on all edges.
// Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
// PDE is treated with central differences on a uniform M x M grid.
// The values of u at the interior points satisfy daes, and
// equations u = 0 at the boundaries are appended, to form a DAE
// system of size N = M^2. Here M = 10.
//
// The system is solved with IDA using the Krylov linear solver
// IDASPGMR. The preconditioner uses the diagonal elements of the
// Jacobian only. Routines for preconditioning, required by
// IDASPGMR, are supplied here. The constraints u >= 0 are posed
// for all components. Output is taken at t = 0, .01, .02, .04,
// ..., 10.24. Two cases are run -- with the Gram-Schmidt type
// being Modified in the first case, and Classical in the second.
// The second run uses IDAReInit and IDAReInitSpgmr.

#include "core/polymec.h"
#include "integrators/dae_integrator.h"
#include "integrators/curtis_powell_reed_preconditioners.h"

/* Problem Constants */

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID
#define ZERO  0.0
#define ONE   1.0
#define TWO   2.0
#define FOUR  4.0

typedef struct 
{  
  long int mm;  // number of grid points 
  real_t dx;
  real_t coeff;

  // Sparsity graph.
  adj_graph_t* sparsity;
} heat2d_t;

// Newly initialized data context.
heat2d_t* heat2d_new()
{
  heat2d_t* data = malloc(sizeof(heat2d_t));
  data->mm  = MGRID;
  data->dx = ONE/(MGRID-ONE);
  data->coeff = ONE/(data->dx * data->dx);

  // Construct a sparsity graph.
  data->sparsity = adj_graph_new(MPI_COMM_SELF, MGRID*MGRID);
  for (long int j = 1; j < MGRID-1; j++) 
  {
    long int offset = data->mm*j;
    for (long int i = 1; i < data->mm-1; i++) 
    {
      long int loc = offset + i;
      int edges[4], num_edges = 4;
      edges[0] = loc-1;
      edges[1] = loc+1;
      edges[2] = loc-data->mm;
      edges[3] = loc+data->mm;

      // Set the edges within the sparsity graph.
      adj_graph_set_num_edges(data->sparsity, loc, num_edges);
      memcpy(adj_graph_edges(data->sparsity, loc), edges, sizeof(int) * num_edges);
    }
  }

  return data;
}

static void heat2d_dtor(void* context)
{
  heat2d_t* data = context;
  adj_graph_free(data->sparsity);
  free(data);
}

static int heat2d_res(void* context, real_t t, real_t* u, real_t* u_dot, real_t* F)
{
  heat2d_t* data = context;
  real_t coeff = data->coeff;
  long int mm = data->mm;
  
  // Initialize F to u, to take care of boundary equations. 
  for (long int j = 0; j < MGRID*MGRID; ++j) 
    F[j] = u[j];
  
  // Loop over interior points; set res = up - (central difference). 
  for (long int j = 1; j < MGRID-1; j++) 
  {
    long int offset = mm*j;
    for (long int i = 1; i < mm-1; i++) 
    {
      long int loc = offset + i;
      real_t dif1 = u[loc-1]  + u[loc+1]  - TWO * u[loc];
      real_t dif2 = u[loc-mm] + u[loc+mm] - TWO * u[loc];
      F[loc]= u_dot[loc] - coeff * (dif1 + dif2);
    }
  }

  return 0;
}

static void heat2d_set_constraints(void* context, real_t* constraints)
{
  // All solution values are non-negative.
  for (int i = 0; i < NEQ; ++i)
    constraints[i] = 1.0;
}

void heat2d_set_initial_conditions(dae_integrator_t* integ, real_t** u, real_t** u_dot)
{
  heat2d_t* data = dae_integrator_context(integ);
  long int mm = data->mm;
  *u = malloc(sizeof(real_t) * NEQ);
  *u_dot = malloc(sizeof(real_t) * NEQ);
  real_t* F = malloc(sizeof(real_t) * NEQ);

  // Initialize u on all grid points. 
  long int mm1 = mm - 1;
  for (long int j = 0; j < mm; j++) 
  {
    real_t yfact = data->dx * j;
    long int offset = mm*j;
    for (long int i = 0;i < mm; i++) 
    {
      real_t xfact = data->dx * i;
      long int loc = offset + i;
      (*u)[loc] = 16.0 * xfact * (ONE - xfact) * yfact * (ONE - yfact);
    }
  }
  
  // Initialize up vector to 0. 
  for (long int i = 0; i < NEQ; ++i)
    (*u_dot)[i] = 0.0;

  // heat2d_res sets F to negative of DAE RHS values at interior points. 
  heat2d_res(data, ZERO, *u, *u_dot, F);

  // Copy -res into u_dot to get correct interior initial up values. 
  for (long int i = 0; i < NEQ; ++i)
    (*u_dot)[i] = -F[i];

  // Set up at boundary points to zero.
  for (long int j = 0; j < mm; j++) 
  {
    long int offset = mm*j;
    for (long int i = 0; i < mm; i++) 
    {
      long int loc = offset + i;
      if (j == 0 || j == mm1 || i == 0 || i == mm1 ) 
        (*u_dot)[loc] = ZERO;
    }
  }

  free(F);
}

// Constructor for heat2d integrator with no preconditioner.
static dae_integrator_t* heat2d_integrator_new()
{
  // Set up a time integrator using GMRES with a Krylov space of maximum 
  // dimension 5.
  heat2d_t* data = heat2d_new();
  dae_integrator_vtable vtable = {.residual = heat2d_res,
                                  .set_constraints = heat2d_set_constraints,
                                  .dtor = heat2d_dtor};
  dae_integrator_t* integ = gmres_dae_integrator_new("heat2d",
                                                     data,
                                                     MPI_COMM_SELF,
                                                     NEQ,
                                                     vtable, 5, 5);

  return integ;
}

// Constructor for block-Jacobi-preconditioned heat2d integrator.
dae_integrator_t* block_jacobi_precond_heat2d_integrator_new()
{
  dae_integrator_t* integ = heat2d_integrator_new();
  heat2d_t* data = dae_integrator_context(integ);
  preconditioner_t* precond = block_jacobi_preconditioner_from_dae_function("Heat 2D", data, heat2d_res, NULL, data->sparsity, NEQ, 1);
  dae_integrator_set_preconditioner(integ, precond);
  return integ;
}

// Constructor for LU-preconditioned heat2d integrator.
dae_integrator_t* lu_precond_heat2d_integrator_new()
{
  dae_integrator_t* integ = heat2d_integrator_new();
  heat2d_t* data = dae_integrator_context(integ);
  preconditioner_t* precond = lu_preconditioner_from_dae_function("Heat 2D", data, heat2d_res, NULL, data->sparsity, NEQ, 1);
  dae_integrator_set_preconditioner(integ, precond);
  return integ;
}

// Constructor for ILU-preconditioned heat2d integrator.
dae_integrator_t* ilu_precond_heat2d_integrator_new()
{
  dae_integrator_t* integ = heat2d_integrator_new();
  ilu_params_t* ilu_params = ilu_params_new();
  heat2d_t* data = dae_integrator_context(integ);
  preconditioner_t* precond = ilu_preconditioner_from_dae_function("Heat 2D", data, heat2d_res, NULL, data->sparsity, NEQ, 1, ilu_params);
  dae_integrator_set_preconditioner(integ, precond);
  return integ;
}

