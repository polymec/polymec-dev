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
  multicomp_poly_basis_t* basis; 
  int basis_dim, num_comp;
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
                               multicomp_poly_basis_t* poly_basis,
                               point_weight_function_t* W)
{
  ASSERT(vtable.num_nodes != NULL);
  ASSERT(vtable.get_nodes != NULL);
  ASSERT(vtable.get_points != NULL);

  gmls_matrix_t* matrix = polymec_malloc(sizeof(gmls_matrix_t));
  matrix->name = string_dup(name);
  matrix->context = context;
  matrix->vtable = vtable;
  matrix->basis = poly_basis;
  matrix->basis_dim = multicomp_poly_basis_dim(poly_basis);
  if (vtable.compute_weight_displacement == NULL)
    matrix->vtable.compute_weight_displacement = simple_weight_displacement;
  matrix->W = W;
  matrix->num_comp = multicomp_poly_basis_num_comp(poly_basis);

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
                               point_t* xjs, int num_nodes, 
                               real_t* W, real_t* Pt, real_t* phi)
{
  int basis_dim = matrix->basis_dim;

  // Compute the moment matrix PtWP.

  // First form Pt*W.
  real_t PtW[basis_dim*num_nodes]; 
  for (int i = 0; i < basis_dim; ++i)
    for (int j = 0; j < num_nodes; ++j)
      PtW[j*basis_dim+i] = Pt[j*basis_dim+i] * W[j]; // PtW_ij = Pt_ij * W_j

  // Now Pt*W*P.
  char no_trans = 'N', trans = 'T';
  real_t alpha = 1.0, beta = 0.0;
  real_t PtWP[basis_dim*basis_dim];
  rgemm(&no_trans, &trans, &basis_dim, &basis_dim, &num_nodes, &alpha, PtW, 
        &basis_dim, Pt, &basis_dim, &beta, PtWP, &basis_dim);
//printf("Pt*W*P = ");
//matrix_fprintf(PtWP, basis_dim, basis_dim, stdout);

  // Now form the matrix phi = (PtWP)^-1*PtW.

  // Factor PtWP all Cholesky-like, since it should be a symmetric matrix.
  char uplo = 'L';
  int info;
  rpotrf(&uplo, &basis_dim, PtWP, &basis_dim, &info); 
  if (info != 0)
  {
    polymec_error("gmls_matrix: Cholesky factorization of P*W*Pt failed. This often means\n"
                  "gmls_matrix: that something is wrong with your point distribution.");
  }

  // Compute (PtWP)^-1 * PtW.
  int nrhs = num_nodes;
  memcpy(phi, PtW, sizeof(real_t) * basis_dim * num_nodes);
  rpotrs(&uplo, &basis_dim, &nrhs, PtWP, &basis_dim, phi, &basis_dim, &info);
}

static void get_neighborhood(gmls_matrix_t* matrix, 
                             int i,
                             point_t* xi, 
                             int* js,
                             point_t* xjs, 
                             int num_nodes,
                             real_t* dx)
{
  int basis_dim = matrix->basis_dim;

  // Get the nodes within this subdomain. Note that we must have 
  // basis_dim <= num_nodes for the matrix to be nonsingular.
  if (num_nodes < basis_dim)
  {
    polymec_error("gmls_matrix: Singular moment matrix for subdomain %d!\n"
                  "gmls_matrix: Number of neighbor nodes N (%d) < polynomial basis dim Q (%d).\n"
                  "gmls_matrix: Nonsingular matrix requires N >= Q.", i, num_nodes, basis_dim);
  }
  ASSERT(num_nodes >= basis_dim); // for nonsingular matrix.
  matrix->vtable.get_nodes(matrix->context, i, js);
  matrix->vtable.get_points(matrix->context, &i, 1, xi);
  matrix->vtable.get_points(matrix->context, js, num_nodes, xjs);

  // Determine the average nodal spacing in the neighborhood.
  *dx = 0;
  int n_dx = 0;
  for (int j = 0; j < num_nodes; ++j)
  {
    real_t D = point_distance(xi, &xjs[j]);
    if ((D == 0.0) && (i != js[j]))
      polymec_error("gmls_matrix: Nodes %d and %d are collocated!", i, js[j]);
    else
    {
      *dx += D;
      ++n_dx;
    }
  }
  ASSERT(n_dx > 0);
  *dx /= n_dx;
}

static void compute_matrix_row(gmls_matrix_t* matrix, 
                               int component,
                               int i,
                               point_t* xi,
                               int* js,
                               point_t* xjs,
                               int num_nodes,
                               real_t* lambdas, 
                               int* columns, 
                               real_t* coeffs)
{
  int basis_dim = matrix->basis_dim;

  // Compute the matrix [Pt]_ij = pi(xj), in column major order.
  real_t Pt[basis_dim*num_nodes];
  for (int j = 0; j < num_nodes; ++j)
  {
    multicomp_poly_basis_compute(matrix->basis, component, 0, 0, 0, 
                                 &xjs[j], &Pt[j*basis_dim]);
  }
//if (i == 0)
//{
//printf("Pt = ");
//matrix_fprintf(Pt, basis_dim, num_nodes, stdout);
//}

  // Compute the (single-component) diagonal matrix W of MLS weights.
  real_t W[num_nodes];
  for (int j = 0; j < num_nodes; ++j)
  {
    vector_t y;
    matrix->vtable.compute_weight_displacement(matrix->context, i, xi, j, &xjs[j], &y);
    W[j] = point_weight_function_value(matrix->W, &y);
  }
//if (i == 0)
//{
//printf("W = ");
//vector_fprintf(W, num_nodes, stdout);
//}

  // Compute the "phi" matrix phi = (Pt * W * P)^-1 * Pt * W.
  real_t phi[basis_dim*num_nodes];
  compute_phi_matrix(matrix, xjs, num_nodes, W, Pt, phi);

  // Now form the coefficients in the matrix row from the product of the 
  // lambda and phi matrices.
  int num_comp = matrix->num_comp;
  memset(coeffs, 0, sizeof(real_t) * num_comp * num_nodes);
  for (int i = 0; i < basis_dim; ++i)
    for (int j = 0; j < num_nodes; ++j)
      coeffs[num_nodes * component + j] += lambdas[i] * phi[basis_dim*j+i];

  // Fill in the column indices.
  for (int n = 0; n < num_nodes; ++n)
    for (int c = 0; c < matrix->num_comp; ++c)
      columns[matrix->num_comp*n+c] = matrix->num_comp * js[n] + c;
}

void gmls_matrix_compute_row(gmls_matrix_t* matrix,
                             int row, 
                             gmls_functional_t* lambda,
                             real_t t,
                             real_t* solution,
                             int* columns,
                             real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  int basis_dim = matrix->basis_dim;
  int component = row % matrix->num_comp;
  int i = row / matrix->num_comp; // index of subdomain

  // Get the neighbor of node i.
  int num_nodes = matrix->vtable.num_nodes(matrix->context, i);
  int js[num_nodes];
  point_t xi, xjs[num_nodes];
  real_t dx;
  get_neighborhood(matrix, i, &xi, js, xjs, num_nodes, &dx);

  // Shift / scale our polynomial basis.
  multicomp_poly_basis_shift(matrix->basis, &xi);
  multicomp_poly_basis_scale(matrix->basis, 1.0/dx);

  // Compute the values of the functional.
  real_t lambdas[matrix->num_comp*basis_dim];
  gmls_functional_compute(lambda, component, i, t, matrix->basis, solution, lambdas);

  // Compute the row using these functional values.
  compute_matrix_row(matrix, component, i, &xi, js, xjs, num_nodes, lambdas, columns, coeffs);
}

void gmls_matrix_compute_dirichlet_row(gmls_matrix_t* matrix,
                                       int row, 
                                       int* columns,
                                       real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  int basis_dim = matrix->basis_dim;
  int component = row % matrix->num_comp;
  int i = row / matrix->num_comp; 

  // Get the neighbor of node i.
  int num_nodes = matrix->vtable.num_nodes(matrix->context, i);
  int js[num_nodes];
  point_t xi, xjs[num_nodes];
  real_t dx;
  get_neighborhood(matrix, i, &xi, js, xjs, num_nodes, &dx);

  // Shift / scale our polynomial basis.
  multicomp_poly_basis_shift(matrix->basis, &xi);
  multicomp_poly_basis_scale(matrix->basis, 1.0/dx);

  // "Compute" the delta function version of the functional.
  real_t delta[matrix->num_comp*basis_dim];
  memset(delta, 0, sizeof(real_t) * matrix->num_comp * basis_dim);
  delta[0] = 1.0; // constant term is unity.
  compute_matrix_row(matrix, component, i, &xi, js, xjs, num_nodes, delta, columns, coeffs);
}

void gmls_matrix_compute_neumann_row(gmls_matrix_t* matrix,
                                     int row,
                                     st_func_t* op,
                                     int* columns,
                                     real_t* coeffs)
{
}

void gmls_matrix_compute_robin_row(gmls_matrix_t* matrix,
                                   int row,
                                   st_func_t* alpha_op,
                                   st_func_t* beta_op,
                                   int* columns,
                                   real_t* coeffs)
{
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

gmls_matrix_t* stencil_based_gmls_matrix_new(multicomp_poly_basis_t* poly_basis,
                                             point_weight_function_t* W,
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
  return gmls_matrix_new("Stencil-based GMLS matrix", sbm, vtable, poly_basis, W);
}

