// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "core/declare_nd_array.h"
#include "meshless/gmls_matrix.h"

struct gmls_matrix_t 
{
  char *name;
  void* context;
  gmls_matrix_vtable vtable;

  point_weight_function_t* W;
  multicomp_poly_basis_t* basis; 
  bool basis_comps_same;
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
  matrix->basis_comps_same = multicomp_poly_basis_components_equal(poly_basis);
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

int gmls_matrix_num_coeffs(gmls_matrix_t* matrix, int i)
{
  return matrix->num_comp * matrix->num_comp * matrix->vtable.num_nodes(matrix->context, i);
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

static void compute_phi_matrix(gmls_matrix_t* matrix, int component,
                               int i, point_t* xi, 
                               point_t* xjs, int num_nodes, 
                               real_t* phi)
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

static void compute_coeffs_for_identical_bases(gmls_matrix_t* matrix, 
                                               int i,
                                               point_t* xi,
                                               int* js,
                                               point_t* xjs,
                                               int num_nodes,
                                               real_t* lambdas, 
                                               int* rows, 
                                               int* columns, 
                                               real_t* coeffs)
{
  int basis_dim = matrix->basis_dim;
  int num_comp = matrix->num_comp;

  // Compute the "phi" matrix phi = (Pt * W * P)^-1 * Pt * W.
  real_t phi[basis_dim*num_nodes];
  compute_phi_matrix(matrix, 0, i, xi, xjs, num_nodes, phi);

  memset(coeffs, 0, sizeof(real_t) * num_comp * num_nodes * num_comp);
  DECLARE_3D_ARRAY(real_t, co, coeffs, num_comp, num_nodes, num_comp);
  int k = 0;
  for (int c = 0; c < num_comp; ++c)
  {
    // Form the coefficients in the matrix from the product of the 
    // lambda and phi matrices. NOTE: phi is stored in column-major order!
    real_t* lam = &(lambdas[c*num_comp*basis_dim]);
    for (int i = 0; i < basis_dim; ++i)
      for (int j = 0; j < num_nodes; ++j)
        co[c][j][c] += lam[i] * phi[basis_dim*j+i];

    // Fill in the row and column indices.
    for (int n = 0; n < num_nodes; ++n)
    {
      for (int cc = 0; cc < matrix->num_comp; ++cc, ++k)
      {
        rows[k] = matrix->num_comp * i + c;
        columns[k] = matrix->num_comp * js[n] + cc;
      }
    }
  }
}

static void compute_coeffs_for_different_bases(gmls_matrix_t* matrix, 
                                               int i,
                                               point_t* xi,
                                               int* js,
                                               point_t* xjs,
                                               int num_nodes,
                                               real_t* lambdas, 
                                               int* rows, 
                                               int* columns, 
                                               real_t* coeffs)
{
  int basis_dim = matrix->basis_dim;
  int num_comp = matrix->num_comp;

  memset(coeffs, 0, sizeof(real_t) * num_comp * num_nodes * num_comp);
  DECLARE_3D_ARRAY(real_t, co, coeffs, num_comp, num_nodes, num_comp);
  int k = 0;
  for (int c = 0; c < num_comp; ++c)
  {
    // Compute the "phi" matrix phi = (Pt * W * P)^-1 * Pt * W.
    real_t phi[basis_dim*num_nodes];
    compute_phi_matrix(matrix, c, i, xi, xjs, num_nodes, phi);

    // Now form the coefficients in the matrix from the product of the 
    // lambda and phi matrices. NOTE: phi is stored in column-major order!
    real_t* lam = &(lambdas[c*num_comp*basis_dim]);
    for (int i = 0; i < basis_dim; ++i)
      for (int j = 0; j < num_nodes; ++j)
        co[c][j][c] += lam[i] * phi[basis_dim*j+i];

    // Fill in the row and column indices.
    for (int n = 0; n < num_nodes; ++n)
    {
      for (int cc = 0; cc < matrix->num_comp; ++cc, ++k)
      {
        rows[k] = matrix->num_comp * i + c;
        columns[k] = matrix->num_comp * js[n] + cc;
      }
    }
  }
}

void gmls_matrix_compute_coeffs(gmls_matrix_t* matrix,
                                int i,
                                gmls_functional_t* lambda,
                                real_t t,
                                real_t* solution,
                                int* rows,
                                int* columns,
                                real_t* coeffs)
{
  // In this function we use the notation in Mirzaei's 2015 paper on 
  // "A new low-cost meshfree method for two and three dimensional 
  //  problems in elasticity."
  int basis_dim = matrix->basis_dim;

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
  real_t lambdas[matrix->num_comp*matrix->num_comp*basis_dim];
  gmls_functional_compute(lambda, i, t, matrix->basis, solution, lambdas);

  // Now compute the matrix coefficients.
  if (matrix->basis_comps_same)
  {
    compute_coeffs_for_identical_bases(matrix, i, &xi, js, xjs, num_nodes,
                                       lambdas, rows, columns, coeffs);
  }
  else
  {
    compute_coeffs_for_different_bases(matrix, i, &xi, js, xjs, num_nodes,
                                       lambdas, rows, columns, coeffs);
  }
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

gmls_functional_t* gmls_matrix_dirichlet_bc_new(gmls_matrix_t* matrix)
{
  return gmls_matrix_robin_bc_new(matrix, NULL, 1.0, 0.0);
}

gmls_functional_t* gmls_matrix_neumann_bc_new(gmls_matrix_t* matrix,
                                              st_func_t* n)
{
  return gmls_matrix_robin_bc_new(matrix, n, 0.0, 1.0);
}

// Collocation volume integral for Robin BC.
static int colocation_num_quad_points(void* context, int i)
{
  return 1;
}

static void colocation_get_quadrature(void* context, int i, point_t* points, real_t* weights)
{
  gmls_matrix_t* matrix = context;
  matrix->vtable.get_points(matrix->context, &i, 1, points);
  weights[0] = 1.0;
}

volume_integral_t* gmls_matrix_bc_quadrature_new(gmls_matrix_t* matrix)
{
  volume_integral_vtable vtable = {.num_quad_points = colocation_num_quad_points,
                                   .get_quadrature = colocation_get_quadrature};
  return volume_integral_new("GMLS BC integral", matrix, vtable);
}

typedef struct
{
  real_t alpha, beta;
  st_func_t* n;
  int num_comp;
} gmls_robin_t;

static void robin_eval_integrands(void* context, real_t t, 
                                  multicomp_poly_basis_t* basis,
                                  point_t* x, vector_t* n, real_t* solution,
                                  real_t* integrands)
{
  int basis_dim = multicomp_poly_basis_dim(basis);
  ASSERT(basis_dim >= 4);

  gmls_robin_t* robin = context;
  int num_comp = robin->num_comp;
  memset(integrands, 0, sizeof(real_t) * num_comp * basis_dim * num_comp);
  DECLARE_3D_ARRAY(real_t, I, integrands, num_comp, basis_dim, num_comp);
  for (int c = 0; c < num_comp; ++c)
  {
    I[c][0][c] = robin->alpha;
    I[c][1][c] = robin->beta;
    I[c][2][c] = robin->beta;
    I[c][3][c] = robin->beta;
  }
}

gmls_functional_t* gmls_matrix_robin_bc_new(gmls_matrix_t* matrix,
                                            st_func_t* n,
                                            real_t alpha,
                                            real_t beta)
{
  ASSERT((alpha != 0.0) || (beta != 0.0));
  ASSERT((n != NULL) || (beta == 0.0));

  int num_comp = multicomp_poly_basis_num_comp(matrix->basis);

#ifndef NDEBUG
  // This functional only works on the standard polynomial basis!
  int degree = multicomp_poly_basis_degree(matrix->basis);
  multicomp_poly_basis_t* std_basis = standard_multicomp_poly_basis_new(num_comp, degree); 
  ASSERT(multicomp_poly_basis_equals(matrix->basis, std_basis));
#endif

  gmls_functional_vtable vtable = {.eval_integrands = robin_eval_integrands,
                                   .dtor = polymec_free};
  gmls_robin_t* robin = polymec_malloc(sizeof(gmls_robin_t));
  robin->alpha = alpha;
  robin->beta = beta;
  robin->n = n;
  robin->num_comp = num_comp;
  volume_integral_t* Qv = gmls_matrix_bc_quadrature_new(matrix);
  char name[1024];
  if (beta == 0.0)
    snprintf(name, 1023, "Dirichlet BC");
  else if (alpha == 0.0)
    snprintf(name, 1023, "Neumann BC");
  else
    snprintf(name, 1023, "Robin BC (alpha = %g, beta = %g)", alpha, beta);
  return volume_gmls_functional_new(name, robin, vtable, num_comp, Qv);
}
