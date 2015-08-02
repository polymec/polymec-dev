// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/gmls_functional.h"

struct gmls_functional_t 
{
  char *name;
  void* context;
  gmls_functional_vtable vtable;

  multicomp_poly_basis_t* basis;
  point_weight_function_t* W;
  int dim, num_comp;
};

gmls_functional_t* gmls_functional_new(const char* name,
                                       void* context,
                                       gmls_functional_vtable vtable,
                                       multicomp_poly_basis_t* poly_basis,
                                       point_weight_function_t* W)
{
  ASSERT(vtable.num_nodes != NULL);
  ASSERT(vtable.get_nodes != NULL);
  ASSERT(vtable.get_points != NULL);
  ASSERT(vtable.num_quad_points != NULL);
  ASSERT(vtable.get_quadrature != NULL);
  ASSERT(vtable.eval_integrands != NULL);

  gmls_functional_t* functional = polymec_malloc(sizeof(gmls_functional_t));
  functional->name = string_dup(name);
  functional->context = context;
  functional->vtable = vtable;
  functional->dim = multicomp_poly_basis_dim(poly_basis);
  functional->basis = poly_basis;
  functional->num_comp = multicomp_poly_basis_num_comp(poly_basis);
  functional->W = W;
  return functional;
}

void gmls_functional_free(gmls_functional_t* functional)
{
  if ((functional->context != NULL) && (functional->vtable.dtor != NULL))
    functional->vtable.dtor(functional->context);
  point_weight_function_free(functional->W);
  polymec_free(functional->name);
  polymec_free(functional);
}

int gmls_functional_num_components(gmls_functional_t* functional)
{
  return functional->num_comp;
}

int gmls_functional_num_nodes(gmls_functional_t* functional, int i)
{
  return functional->vtable.num_nodes(functional->context, i);
}

static void compute_simple_Phi_matrix(gmls_functional_t* functional,
                                      point_t* xjs, int num_nodes,
                                      real_t* W_sc, real_t* Phi)
{
  int num_comp = functional->num_comp;
  int Q = functional->dim;

  // Compute the single­component matrix [P]_ij = pi(xj), in column major order.
  real_t P_sc[num_nodes*Q];
  for (int j = 0; j < num_nodes; ++j)
    multicomp_poly_basis_compute(functional->basis, 0, 0, 0, 0, &xjs[j], &P_sc[j*Q]);

  // Compute the single-component moment matrix PWPt_sc.
  real_t WPt_sc[num_nodes*Q], PWPt_sc[Q*Q];
  char no_trans = 'N', trans = 'T';
  real_t alpha = 1.0, beta = 0.0;
  rgemm(&no_trans, &trans, &num_nodes, &Q, &num_nodes, &alpha, W_sc, 
        &num_nodes, P_sc, &Q, &beta, WPt_sc, &num_nodes);
  rgemm(&no_trans, &no_trans, &Q, &Q, &num_nodes, &alpha, P_sc, 
        &Q, WPt_sc, &num_nodes, &beta, PWPt_sc, &Q);

  // Now form the single-component Phi matrix Phi_sc = (PWPt)^-1 * Pt * W.

  // Factor PWPt.
  char uplo = 'L';
  int info;
  rpotrf(&uplo, &Q, PWPt_sc, &Q, &info); 
  ASSERT(info == 0);

  // Compute (PWPt)^1 * PtW.
  real_t Phi_sc[Q*num_nodes];
  rgemm(&trans, &no_trans, &Q, &num_nodes, &num_nodes, &alpha, P_sc, 
        &num_nodes, W_sc, &num_nodes, &beta, Phi_sc, &Q);
  int one = 1;
  rpotrs(&uplo, &Q, &one, PWPt_sc, &Q, Phi_sc, &Q, &info);

  // Fill in Phi with block diagonal Phi_sc's.
  memset(Phi, 0, sizeof(real_t) * num_comp * Q * num_comp * num_nodes);
  for (int c = 0; c < num_comp; ++c)
    memcpy(&Phi[num_comp*Q+num_comp*c], Phi_sc, sizeof(real_t) * Q * num_nodes);
}

static void compute_general_Phi_matrix(gmls_functional_t* functional,
                                       point_t* xjs, int num_nodes,
                                       real_t* W_sc, real_t* Phi)
{
  // Compute the matrix [P]_ij = pi(xj), in column major order.
  // FIXME: this has richer structure in general than in Mirzaei's 2015 paper
  // FIXME: on elasticity!
  POLYMEC_NOT_IMPLEMENTED;
}

void gmls_functional_compute_block_row(gmls_functional_t* functional,
                                       real_t t,
                                       int i,
                                       int* rows,
                                       int* columns,
                                       real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  int num_comp = functional->num_comp;

  // Compute the quadrature points/weights for this subdomain.
  int num_quad_points = functional->vtable.num_quad_points(functional->context, i);
  point_t quad_points[num_quad_points];
  real_t quad_weights[num_quad_points];
  functional->vtable.get_quadrature(functional->context, i, num_quad_points, 
                             quad_points, quad_weights);

  // Loop through the points and compute the lambda matrix of functional 
  // approximants.
  int Q = functional->dim;
  real_t lambda[num_comp*num_comp*Q];
  memset(lambda, 0, sizeof(real_t) * num_comp*num_comp*Q);
  for (int q = 0; q < num_quad_points; ++q)
  {
    point_t* xq = &quad_points[q];
    real_t wq = quad_weights[q];

    // Now compute the (multi-component) integrands for the functional at 
    // this point.
    real_t integrands[num_comp*num_comp*Q];
    functional->vtable.eval_integrands(functional->context, t, xq, functional->basis, integrands);

    // Integrate.
    for (int i = 0; i < num_comp*num_comp*Q; ++i)
      lambda[i] += wq * integrands[i];
  }

  // Get the nodes within this subdomain.
  int num_nodes = functional->vtable.num_nodes(functional->context, i);
  int nodes[num_nodes];
  functional->vtable.get_nodes(functional->context, i, nodes);
  point_t xi, xjs[num_nodes];
  functional->vtable.get_points(functional->context, &i, 1, &xi);
  functional->vtable.get_points(functional->context, nodes, num_nodes, &xi);

  // Compute the single-component diagonal matrix of MLS weights.
  real_t W_sc[num_nodes*num_nodes];
  memset(W_sc, 0.0, sizeof(real_t) * num_nodes*num_nodes);
  for (int j = 0; j < num_nodes; ++j)
  {
    vector_t y;
    point_displacement(&xi, &xjs[j], &y);
    W_sc[num_nodes*j+j] = point_weight_function_value(functional->W, &y);
  }

  // Compute the "Phi" matrix Phi = (Pt * W * P)^-1 * W * Pt. (This is simpler 
  // if all the components of the polynomial basis are equal.)
  real_t Phi[num_comp*Q*num_comp*num_nodes];
  if (multicomp_poly_basis_components_equal(functional->basis))
    compute_simple_Phi_matrix(functional, xjs, num_nodes, W_sc, Phi);
  else
    compute_general_Phi_matrix(functional, xjs, num_nodes, W_sc, Phi);

  // Now form the matrix coefficients from the product of the lambda and Phi
  // matrices.
  char no_trans = 'N';
  int M = num_comp; // Rows in lambda matrix.
  int N = num_comp * num_nodes; // Columns in Phi matrix.
  int K = num_comp * Q; // Columns in lambda, rows in Phi.
  real_t alpha = 1.0, beta = 0.0;
  rgemm(&no_trans, &no_trans, &M, &N, &K, &alpha, lambda, &M, Phi, &K, 
        &beta,coeffs, &M);

  // Fill in the row and column indices.
  for (int c1 = 0; c1 < num_comp; ++c1)
  {
    for (int n = 0; n < num_nodes; ++n)
    {
      for (int c2 = 0; c2 < num_comp; ++c2)
      {
        int j = c1*num_comp*num_nodes + num_comp*n + c2;
        rows[j] = num_comp * i + c1;
        columns[j] = num_comp * nodes[n] + c2;
      }
    }
  }
}

