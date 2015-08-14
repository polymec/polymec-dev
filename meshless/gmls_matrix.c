// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "meshless/gmls_matrix.h"

struct gmls_matrix_t 
{
  char *name;
  void* context;
  gmls_matrix_vtable vtable;

  point_weight_function_t* W;
  int num_comp;
};

static void simple_weight_displacement(void* context, 
                                       int i,
                                       point_t* xi, 
                                       int j,
                                       point_t* xj,
                                       vector_t* y)
{
  point_displacement(xi, xj, y);
}
  
gmls_matrix_t* gmls_matrix_new(const char* name,
                               void* context,
                               gmls_matrix_vtable vtable,
                               point_weight_function_t* W,
                               int num_components)
{
  ASSERT(vtable.num_nodes != NULL);
  ASSERT(vtable.get_nodes != NULL);
  ASSERT(vtable.get_points != NULL);
  ASSERT(num_components > 0);

  gmls_matrix_t* matrix = polymec_malloc(sizeof(gmls_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  if (vtable.compute_weight_displacement == NULL)
    matrix->vtable.compute_weight_displacement = simple_weight_displacement;
  matrix->W = W;
  matrix->num_comp = num_components;

  return matrix;
}

void gmls_matrix_free(gmls_matrix_t* matrix)
{
  if ((matrix->context != NULL) && (matrix->vtable.dtor != NULL))
    matrix->vtable.dtor(matrix->context);
  point_weight_function_free(matrix->W);
  polymec_free(matrix->name);
  polymec_free(matrix);
}

int gmls_matrix_num_columns(gmls_matrix_t* matrix, int row)
{
  int i = row / matrix->num_comp; // index of subdomain
  return matrix->num_comp * matrix->vtable.num_nodes(matrix->context, i);
}

static void compute_phi_matrix(gmls_matrix_t* matrix,
                               point_t* xjs, int num_nodes, int basis_dim,
                               real_t* W, real_t* P, real_t* phi)
{
  // Compute the moment matrix PWPt.
  real_t WPt[num_nodes*basis_dim], PWPt[basis_dim*basis_dim];
  char no_trans = 'N', trans = 'T';
  real_t alpha = 1.0, beta = 0.0;
  rgemm(&no_trans, &trans, &num_nodes, &basis_dim, &num_nodes, &alpha, W, 
        &num_nodes, P, &basis_dim, &beta, WPt, &num_nodes);
  rgemm(&no_trans, &no_trans, &basis_dim, &basis_dim, &num_nodes, &alpha, P, 
        &basis_dim, WPt, &num_nodes, &beta, PWPt, &basis_dim);
matrix_fprintf(P, basis_dim, num_nodes, stdout);
printf("\n\n\n");
matrix_fprintf(W, num_nodes, num_nodes, stdout);
printf("\n\n\n");
matrix_fprintf(PWPt, basis_dim, basis_dim, stdout);

  // Now form the matrix phi = (PWPt)^-1 * Pt * W.

  // Factor PWPt all Cholesky-like, since it should be a symmetric matrix.
  char uplo = 'L';
  int info;
  rpotrf(&uplo, &basis_dim, PWPt, &basis_dim, &info); 
  ASSERT(info == 0);

  // Compute (PWPt)^1 * PtW.
  rgemm(&trans, &no_trans, &basis_dim, &num_nodes, &num_nodes, &alpha, P, 
        &num_nodes, W, &num_nodes, &beta, phi, &basis_dim);
  int one = 1;
  rpotrs(&uplo, &basis_dim, &one, PWPt, &basis_dim, phi, &basis_dim, &info);
}

static void compute_matrix_row(gmls_matrix_t* matrix, 
                               multicomp_poly_basis_t* poly_basis,
                               int component,
                               int i, 
                               real_t* lambdas, 
                               int* columns, 
                               real_t* coeffs)
{
  int basis_dim = multicomp_poly_basis_dim(poly_basis);

  // Get the nodes within this subdomain. Note that we must have 
  // basis_dim <= num_nodes for the matrix to be nonsingular.
  int num_nodes = matrix->vtable.num_nodes(matrix->context, i);

  if (num_nodes < basis_dim)
    polymec_error("gmls_matrix: Singular moment matrix!\n"
                  "gmls_matrix: Number of neighbor nodes N (%d) < polynomial basis dim Q (%d).\n"
                  "gmls_matrix: Nonsingular matrix requires N >= Q.", num_nodes, basis_dim);
  ASSERT(num_nodes >= basis_dim); // for nonsingular matrix.
  int nodes[num_nodes];
  matrix->vtable.get_nodes(matrix->context, i, nodes);
  point_t xi, xjs[num_nodes];
  matrix->vtable.get_points(matrix->context, &i, 1, &xi);
  matrix->vtable.get_points(matrix->context, nodes, num_nodes, xjs);

  // Compute the matrix [P]_ij = pi(xj), in column major order.
  real_t P[basis_dim*num_nodes];
  for (int j = 0; j < num_nodes; ++j)
  {
    multicomp_poly_basis_compute(poly_basis, component, 0, 0, 0, 
                                 &xjs[j], &P[j*basis_dim]);
  }

  // Compute the (single-component) diagonal matrix W of MLS weights.
  real_t W[num_nodes*num_nodes];
  memset(W, 0.0, sizeof(real_t) * num_nodes*num_nodes);
  for (int j = 0; j < num_nodes; ++j)
  {
    vector_t y;
    matrix->vtable.compute_weight_displacement(matrix->context, i, &xi, j, &xjs[j], &y);
    W[num_nodes*j+j] = point_weight_function_value(matrix->W, &y);
  }

  // Compute the "Phi" matrix Phi = (Pt * W * P)^-1 * W * Pt.
  real_t phi[basis_dim*num_nodes];
  compute_phi_matrix(matrix, xjs, num_nodes, basis_dim, W, P, phi);

  // Now form the matrix coefficients from the product of the lambda and Phi
  // matrices.
  char no_trans = 'N';
  int M = 1; // Rows in lambda matrix.
  int N = num_nodes; // Columns in Phi matrix.
  int K = matrix->num_comp*basis_dim; // Columns in lambda, rows in Phi.
  real_t alpha = 1.0, beta = 0.0;
  rgemm(&no_trans, &no_trans, &M, &N, &K, &alpha, lambdas, &M, phi, &K, 
        &beta, coeffs, &M);

  // Fill in the column indices.
  for (int n = 0; n < num_nodes; ++n)
    for (int c = 0; c < matrix->num_comp; ++c)
      columns[matrix->num_comp*n+c] = matrix->num_comp * nodes[n] + c;
}

void gmls_matrix_compute_row(gmls_matrix_t* matrix,
                             int row, 
                             gmls_functional_t* lambda,
                             real_t t,
                             int* columns,
                             real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  multicomp_poly_basis_t* poly_basis = gmls_functional_basis(lambda);
  int basis_dim = multicomp_poly_basis_dim(poly_basis);
  ASSERT(matrix->num_comp == multicomp_poly_basis_num_comp(poly_basis));

  // Compute the values of the functional.
  int component = row % matrix->num_comp;
  int i = row / matrix->num_comp; // index of subdomain
  real_t lambdas[matrix->num_comp*basis_dim];
  gmls_functional_compute(lambda, component, i, t, lambdas);

  // Compute the row using these functional values.
  compute_matrix_row(matrix, poly_basis, component, i, lambdas, columns, coeffs);
}

void gmls_matrix_compute_dirichlet_row(gmls_matrix_t* matrix,
                                       int row, 
                                       gmls_functional_t* lambda,
                                       int* columns,
                                       real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  multicomp_poly_basis_t* poly_basis = gmls_functional_basis(lambda);
  int basis_dim = multicomp_poly_basis_dim(poly_basis);
  ASSERT(matrix->num_comp == multicomp_poly_basis_num_comp(poly_basis));

  int component = row % matrix->num_comp;
  int i = row / matrix->num_comp; 
  real_t delta[matrix->num_comp*basis_dim];
  memset(delta, 0, sizeof(real_t) * matrix->num_comp * basis_dim);
  for (int j = 0; j < basis_dim; ++j)
    delta[matrix->num_comp*j+component] = 1.0;
  compute_matrix_row(matrix, poly_basis, component, i, delta, columns, coeffs);
}

// Stencil-based GMLS matrix.
typedef struct
{
  point_cloud_t* points;
  real_t* extents;
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

static void sbm_compute_y(void* context, 
                          int i, point_t* xi, 
                          int j, point_t* xj, 
                          vector_t* y)
{
  sbm_t* sbm = context;
  real_t Ri = sbm->extents[i];
  y->x = (xj->x - xi->x) / Ri;
  y->y = (xj->y - xi->y) / Ri;
  y->z = (xj->z - xi->z) / Ri;
}

static void sbm_dtor(void* context)
{
  sbm_t* sbm = context;
  polymec_free(sbm);
}

gmls_matrix_t* stencil_based_gmls_matrix_new(point_weight_function_t* W,
                                             int num_components,
                                             point_cloud_t* points,
                                             real_t* extents,
                                             stencil_t* stencil)
{
  sbm_t* sbm = polymec_malloc(sizeof(sbm_t));
  sbm->points = points;
  sbm->extents = extents;
  sbm->stencil = stencil;
  gmls_matrix_vtable vtable = {.num_nodes = sbm_num_nodes,
                               .get_nodes = sbm_get_nodes,
                               .get_points = sbm_get_points,
                               .compute_weight_displacement = sbm_compute_y,
                               .dtor = sbm_dtor};
  return gmls_matrix_new("Stencil-based GMLS matrix", sbm, vtable, W, num_components);
}

