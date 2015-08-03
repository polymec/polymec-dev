// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/gmls_matrix.h"

struct gmls_matrix_t 
{
  char *name;
  void* context;
  gmls_matrix_vtable vtable;

  gmls_functional_t* lambda;
  point_weight_function_t* W;

  multicomp_poly_basis_t* basis;
  int dim, num_comp;
};

gmls_matrix_t* gmls_matrix_new(const char* name,
                               void* context,
                               gmls_matrix_vtable vtable,
                               gmls_functional_t* lambda,
                               point_weight_function_t* W)
{
  ASSERT(vtable.num_nodes != NULL);
  ASSERT(vtable.get_nodes != NULL);
  ASSERT(vtable.get_points != NULL);

  gmls_matrix_t* matrix = polymec_malloc(sizeof(gmls_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  matrix->lambda = lambda;
  matrix->W = W;

  // Jot some stuff down for quick reference.
  matrix->basis = gmls_functional_basis(lambda);
  matrix->dim = multicomp_poly_basis_dim(matrix->basis);
  matrix->num_comp = multicomp_poly_basis_num_comp(matrix->basis);
  return matrix;
}

void gmls_matrix_free(gmls_matrix_t* matrix)
{
  if ((matrix->context != NULL) && (matrix->vtable.dtor != NULL))
    matrix->vtable.dtor(matrix->context);
  point_weight_function_free(matrix->W);
  gmls_functional_free(matrix->lambda);
  matrix->basis = NULL;
  polymec_free(matrix->name);
  polymec_free(matrix);
}

gmls_functional_t* gmls_matrix_functional(gmls_matrix_t* matrix)
{
  return matrix->lambda;
}

int gmls_matrix_num_columns(gmls_matrix_t* matrix, int row)
{
  int i = row / matrix->num_comp; // index of subdomain
  return matrix->num_comp * matrix->vtable.num_nodes(matrix->context, i);
}

static void compute_phi_matrix(gmls_matrix_t* matrix,
                               point_t* xjs, int num_nodes,
                               real_t* W, real_t* phi)
{
  int Q = matrix->dim;

  // Compute the matrix [P]_ij = pi(xj), in column major order.
  real_t P[num_nodes*Q];
  for (int j = 0; j < num_nodes; ++j)
    multicomp_poly_basis_compute(matrix->basis, 0, 0, 0, 0, &xjs[j], &P[j*Q]);

  // Compute the moment matrix PWPt.
  real_t WPt[num_nodes*Q], PWPt[Q*Q];
  char no_trans = 'N', trans = 'T';
  real_t alpha = 1.0, beta = 0.0;
  rgemm(&no_trans, &trans, &num_nodes, &Q, &num_nodes, &alpha, W, 
        &num_nodes, P, &Q, &beta, WPt, &num_nodes);
  rgemm(&no_trans, &no_trans, &Q, &Q, &num_nodes, &alpha, P, 
        &Q, WPt, &num_nodes, &beta, PWPt, &Q);

  // Now form the matrix phi = (PWPt)^-1 * Pt * W.

  // Factor PWPt.
  char uplo = 'L';
  int info;
  rpotrf(&uplo, &Q, PWPt, &Q, &info); 
  ASSERT(info == 0);

  // Compute (PWPt)^1 * PtW.
  rgemm(&trans, &no_trans, &Q, &num_nodes, &num_nodes, &alpha, P, 
        &num_nodes, W, &num_nodes, &beta, phi, &Q);
  int one = 1;
  rpotrs(&uplo, &Q, &one, PWPt, &Q, phi, &Q, &info);
}

void gmls_matrix_compute_row(gmls_matrix_t* matrix,
                             int row, 
                             real_t t,
                             int* columns,
                             real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  int num_comp = matrix->num_comp;
  int Q = matrix->dim;

  // Compute the values of the functional.
  int component = row % num_comp;
  int i = row / num_comp; // index of subdomain
  real_t lambdas[Q];
  gmls_functional_compute(matrix->lambda, component, i, t, lambdas);

  // Get the nodes within this subdomain.
  int num_nodes = matrix->vtable.num_nodes(matrix->context, i);
  int nodes[num_nodes];
  matrix->vtable.get_nodes(matrix->context, i, nodes);
  point_t xi, xjs[num_nodes];
  matrix->vtable.get_points(matrix->context, &i, 1, &xi);
  matrix->vtable.get_points(matrix->context, nodes, num_nodes, &xi);

  // Compute the (single-component) diagonal matrix W of MLS weights.
  real_t W[num_nodes*num_nodes];
  memset(W, 0.0, sizeof(real_t) * num_nodes*num_nodes);
  for (int j = 0; j < num_nodes; ++j)
  {
    vector_t y;
    point_displacement(&xi, &xjs[j], &y);
    W[num_nodes*j+j] = point_weight_function_value(matrix->W, &y);
  }

  // Compute the "Phi" matrix Phi = (Pt * W * P)^-1 * W * Pt.
  real_t phi[Q*num_nodes];
  compute_phi_matrix(matrix, xjs, num_nodes, W, phi);

  // Now form the matrix coefficients from the product of the lambda and Phi
  // matrices.
  char no_trans = 'N';
  int M = 1; // Rows in lambda matrix.
  int N = num_nodes; // Columns in Phi matrix.
  int K = Q; // Columns in lambda, rows in Phi.
  real_t alpha = 1.0, beta = 0.0;
  rgemm(&no_trans, &no_trans, &M, &N, &K, &alpha, lambdas, &M, phi, &K, 
        &beta, coeffs, &M);

  // Fill in the column indices.
  for (int n = 0; n < num_nodes; ++n)
    for (int c = 0; c < num_comp; ++c)
      columns[n] = num_comp * nodes[n] + c;
}

// Stencil-based GMLS matrix.
typedef struct
{
  point_cloud_t* points;
  stencil_t* stencil;
} sbm_t;

static int sbm_num_nodes(void* context, int i)
{
  sbm_t* sbm = context;
  return stencil_size(sbm->stencil, i);
}

static void sbm_get_nodes(void* context, int i, int* nodes)
{
  sbm_t* sbm = context;
  stencil_get_neighbors(sbm->stencil, i, nodes);
}

static void sbm_get_points(void* context, int* nodes, int num_nodes, point_t* points)
{
  sbm_t* sbm = context;
  for (int i = 0; i < num_nodes; ++i)
    points[i] = sbm->points->points[nodes[i]];
}

static void sbm_dtor(void* context)
{
  sbm_t* sbm = context;
  polymec_free(sbm);
}

gmls_matrix_t* stencil_based_gmls_matrix_new(gmls_functional_t* lambda,
                                             point_weight_function_t* W,
                                             point_cloud_t* points,
                                             stencil_t* stencil)
{
  sbm_t* sbm = polymec_malloc(sizeof(sbm_t));
  sbm->points = points;
  sbm->stencil = stencil;
  gmls_matrix_vtable vtable = {.num_nodes = sbm_num_nodes,
                               .get_nodes = sbm_get_nodes,
                               .get_points = sbm_get_points,
                               .dtor = sbm_dtor};
  return gmls_matrix_new("Stencil-based GMLS matrix", sbm, vtable, lambda, W);
}

