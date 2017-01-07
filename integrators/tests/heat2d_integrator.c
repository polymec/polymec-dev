// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
#include "core/declare_nd_array.h"
#include "integrators/dae_integrator.h"
#include "integrators/bj_newton_pc.h"

/* Problem Constants */

#define NOUT  11
#define MGRID 10
#define NEQ   MGRID*MGRID

typedef struct 
{  
  long int mm;  // number of grid points 
  real_t dx;
  real_t coeff;

  // Sparsity graph.
  adj_graph_t* sparsity;
} heat2d_t;

// Newly initialized data context.
heat2d_t* heat2d_new(void);
heat2d_t* heat2d_new()
{
  heat2d_t* data = polymec_malloc(sizeof(heat2d_t));
  data->mm  = MGRID;
  data->dx = 1.0/(MGRID-1);
  data->coeff = 1.0/(data->dx * data->dx);

  // Construct a sparsity graph.
  data->sparsity = adj_graph_new(MPI_COMM_SELF, MGRID*MGRID);
  for (int i = 1; i < MGRID-1; ++i) 
  {
    for (int j = 1; j < MGRID-1; ++j) 
    {
      int i_self = ARRAY_INDEX_2D(MGRID, MGRID, i, j);
      int edges[4], num_edges = 4;
      edges[0] = ARRAY_INDEX_2D(MGRID, MGRID, i-1, j);
      edges[1] = ARRAY_INDEX_2D(MGRID, MGRID, i+1, j);
      edges[2] = ARRAY_INDEX_2D(MGRID, MGRID, i, j-1);
      edges[3] = ARRAY_INDEX_2D(MGRID, MGRID, i, j+1);

      // Set the edges within the sparsity graph.
      adj_graph_set_num_edges(data->sparsity, i_self, num_edges);
      memcpy(adj_graph_edges(data->sparsity, i_self), edges, sizeof(int) * num_edges);
    }
  }

  return data;
}

static void heat2d_dtor(void* context)
{
  heat2d_t* data = context;
  adj_graph_free(data->sparsity);
  polymec_free(data);
}

static int heat2d_F(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F)
{
  heat2d_t* data = context;
  real_t coeff = data->coeff;
  
  // Initialize F to u, to take care of boundary equations. 
  for (long int j = 0; j < MGRID*MGRID; ++j) 
    F[j] = U[j];
  
  // Loop over interior points; set res = up - (central difference). 
  DECLARE_2D_ARRAY(real_t, U_ij, U, MGRID, MGRID);
  DECLARE_2D_ARRAY(real_t, U_dot_ij, U_dot, MGRID, MGRID);
  DECLARE_2D_ARRAY(real_t, F_ij, F, MGRID, MGRID);
  for (int i = 1; i < MGRID-1; ++i) 
  {
    for (int j = 1; j < MGRID-1; ++j) 
    {
      real_t dif1 = U_ij[i-1][j] + U_ij[i+1][j] - 2.0*U_ij[i][j];
      real_t dif2 = U_ij[i][j-1] + U_ij[i][j+1] - 2.0*U_ij[i][j];
      F_ij[i][j] = U_dot_ij[i][j] - coeff * (dif1 + dif2);
    }
  }

  return 0;
}

static int heat2d_J(void* context, 
                    real_t t, real_t* U, 
                    real_t alpha, real_t* U_dot, 
                    real_t* F, krylov_matrix_t* J)
{
  heat2d_t* data = context;
  real_t coeff = data->coeff;
  
  // Boundary points correspond to rows of the identity matrix.
  for (int i = 0; i < MGRID; ++i) 
  {
    index_t num_cols = 1, row, indices[1];
    real_t values[1] = {1.0};
    row = indices[0] = ARRAY_INDEX_2D(MGRID, MGRID, 0, i);
    krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
    row = indices[0] = ARRAY_INDEX_2D(MGRID, MGRID, MGRID-1, i);
    krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
    row = indices[0] = ARRAY_INDEX_2D(MGRID, MGRID, i, 0);
    krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
    row = indices[0] = ARRAY_INDEX_2D(MGRID, MGRID, i, MGRID-1);
    krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
  }

  // Loop over all interior grid points. 
  for (int i = 1; i < MGRID-1; ++i) 
  {
    for (int j = 1; j < MGRID-1; ++j) 
    {
      real_t J_self = 0.0, J_left = 0.0, J_right = 0.0, J_up = 0.0, J_down = 0.0;
      index_t I_self  = ARRAY_INDEX_2D(MGRID, MGRID, i, j), 
              I_left  = ARRAY_INDEX_2D(MGRID, MGRID, i-1, j),
              I_right = ARRAY_INDEX_2D(MGRID, MGRID, i+1, j),
              I_up    = ARRAY_INDEX_2D(MGRID, MGRID, i, j+1),
              I_down  = ARRAY_INDEX_2D(MGRID, MGRID, i, j-1);

      // Add in terms. 
      J_self = alpha + 4.0 * coeff;
      J_left = -coeff;
      J_right = -coeff;
      J_up = -coeff;
      J_down = -coeff;

      // Stick 'em into the matrix.
      index_t row = I_self, num_cols = 5;
      index_t indices[5] = {I_self, I_left, I_right, I_up, I_down};
      real_t values[5] = {J_self, J_left, J_right, J_up, J_down};
      krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
    }
  }

  // Assemble the Jacobian matrix.
  krylov_matrix_assemble(J);

  return 0;
}

void heat2d_set_initial_conditions(dae_integrator_t* integ, real_t** U, real_t** U_dot);
void heat2d_set_initial_conditions(dae_integrator_t* integ, real_t** U, real_t** U_dot)
{
  heat2d_t* data;
  const char* integ_name = dae_integrator_name(integ);
  if (string_contains(integ_name, "INK"))
    data = ink_dae_integrator_context(integ);
  else
    data = dae_integrator_context(integ);
  *U = polymec_malloc(sizeof(real_t) * NEQ);
  *U_dot = polymec_malloc(sizeof(real_t) * NEQ);
  real_t* F = polymec_malloc(sizeof(real_t) * NEQ);

  DECLARE_2D_ARRAY(real_t, U_ij, *U, MGRID, MGRID);
  DECLARE_2D_ARRAY(real_t, U_dot_ij, *U_dot, MGRID, MGRID);

  // Initialize u on all grid points. 
  for (int i = 0; i < MGRID; ++i)
  {
    real_t xfact = data->dx * i;
    for (int j = 0; j < MGRID; ++j) 
    {
      real_t yfact = data->dx * j;
      U_ij[i][j] = 16.0 * xfact * (1.0 - xfact) * yfact * (1.0 - yfact);
    }
  }
  
  // Initialize U_dot vector to 0. 
  memset(*U_dot, 0, sizeof(real_t) * NEQ);

  // heat2d_F sets F to negative of DAE RHS values at interior points. 
  heat2d_F(data, 0.0, *U, *U_dot, F);

  // Copy -res into u_dot to get correct interior initial up values. 
  for (int i = 0; i < NEQ; ++i)
    (*U_dot)[i] = -F[i];

  // Set U_dot at boundary points to zero.
  for (int i = 0; i < MGRID; ++i) 
  {
    U_dot_ij[i][0] = 0.0;
    U_dot_ij[i][MGRID-1] = 0.0;
    U_dot_ij[0][i] = 0.0;
    U_dot_ij[MGRID-1][i] = 0.0;
  }

  polymec_free(F);
}

// Constructor for Jacobian-Free Newton-Krylov heat2d integrator with 
// the given preconditioner.
static dae_integrator_t* jfnk_heat2d_integrator_new(heat2d_t* data, newton_pc_t* precond)
{
  // Set up a time integrator using GMRES with a Krylov space of maximum 
  // dimension 5.
  dae_integrator_t* integ = jfnk_dae_integrator_new(5, MPI_COMM_SELF, 
                                                    DAE_ALL_DIFFERENTIAL, 
                                                    DAE_ALL_NONNEGATIVE,
                                                    NEQ, 0, data, 
                                                    heat2d_F, NULL, heat2d_dtor,
                                                    precond, JFNK_DAE_GMRES, 5);

  return integ;
}

// Constructor for block-Jacobi-preconditioned heat2d integrator.
dae_integrator_t* bj_pc_jfnk_heat2d_integrator_new(void);
dae_integrator_t* bj_pc_jfnk_heat2d_integrator_new()
{
  heat2d_t* data = heat2d_new();
  newton_pc_t* precond = dae_cpr_bj_newton_pc_new(MPI_COMM_WORLD, data, heat2d_F, NULL, data->sparsity, NEQ, 0, 1);
  return jfnk_heat2d_integrator_new(data, precond);
}

// Constructor for Inexact Newton-Krylov heat2d integrator.
dae_integrator_t* ink_heat2d_integrator_new(krylov_factory_t* factory);
dae_integrator_t* ink_heat2d_integrator_new(krylov_factory_t* factory)
{
  heat2d_t* data = heat2d_new();
  matrix_sparsity_t* J_sparsity = matrix_sparsity_from_graph(data->sparsity, NULL);

  // Set up a nonlinear solver using GMRES with a full Newton step.
  return ink_dae_integrator_new(5, MPI_COMM_SELF, 
                                DAE_ALL_DIFFERENTIAL, DAE_ALL_NONNEGATIVE,
                                factory, J_sparsity, data,
                                heat2d_F, heat2d_J, heat2d_dtor);
}

