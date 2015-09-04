// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/declare_nd_array.h"
#include "core/linear_algebra.h"
#include "meshless/mlpg_quadrature.h"
#include "elastic_gmls_functional.h"

typedef struct
{
  real_t E, nu;
  int degree;
  point_cloud_t* points;
  real_t* subdomain_extents;
} elastic_t;

// This helper evaluates the coefficients of the stress-strain matrix D, 
// storing them in column-major order in the array D.
static void eval_stress_strain_matrix(real_t E, real_t nu, real_t* D)
{
  real_t k1 = E / ((1.0 - 2.0*nu) * (1.0 + nu));
  real_t k2 = E / (2.0 * (1.0 + nu));
  memset(D, 0, sizeof(real_t) * 36);
  D[0] = k1 * (1.0 - nu); D[6] = k1 * nu        ; D[12] = k1 * nu;
  D[1] = k1 * nu        ; D[7] = k1 * (1.0 - nu); D[13] = k1 * nu;
  D[2] = k1 * nu        ; D[8] = k1 * nu        ; D[14] = k1 * nu;
  D[21] = k2 * 1.0      ; D[28] = k2 * 1.0      ; D[35] = k2 * 1.0;
}

// This helper evaluates the coefficients of the matrix of normal vectors N, 
// storing them in column-major order in the array N.
static void eval_normal_vector_matrix(vector_t* n, real_t* N)
{
  real_t n1 = n->x;
  real_t n2 = n->y;
  real_t n3 = n->z;
  memset(N, 0, sizeof(real_t) * 18);
  N[0] = n1; N[4] = n2; N[8] = n3; 
  N[10] = n3; N[11] = n2; N[12] = n3;
  N[14] = n1; N[15] = n2; N[16] = n1;
}

// This helper evaluates the coefficients of the polynomial basis matrix P, 
// storing them in column-major order in the array P.
static void eval_basis_matrix(real_t dpdx, real_t dpdy, real_t dpdz, real_t* P)
{
  memset(P, 0, sizeof(real_t) * 18);
  P[0] = dpdx;
  P[4] = dpdz;
  P[5] = dpdy;
  P[7] = dpdy;
  P[9] = dpdz;
  P[11] = dpdx;
  P[14] = dpdz;
  P[15] = dpdy;
  P[16] = dpdx;
}

static void elastic_eval_integrands(void* context, real_t t, 
                                    multicomp_poly_basis_t* basis, 
                                    point_t* x, vector_t* n, real_t* solution,
                                    real_t* integrands)
{
  elastic_t* elastic = context;
  int dim = multicomp_poly_basis_dim(basis);

  // Evaluate the polynomial derivatives. All the components use the 
  // same basis, so we don't have to distinguish between them
  real_t dpdx[dim], dpdy[dim], dpdz[dim];
  multicomp_poly_basis_compute(basis, 0, 1, 0, 0, x, dpdx);
  multicomp_poly_basis_compute(basis, 0, 0, 1, 0, x, dpdy);
  multicomp_poly_basis_compute(basis, 0, 0, 0, 1, x, dpdz);

  // Evaluate the D and N matrices.
  real_t D[6*6], N[3*6];
  eval_stress_strain_matrix(elastic->E, elastic->nu, D);
  eval_normal_vector_matrix(n, N);

  memset(integrands, 0, sizeof(int) * dim);
  DECLARE_3D_ARRAY(real_t, I, integrands, 3, 3, dim);
  for (int k = 0; k < dim; ++k)
  {
    // Evaluate the P matrix for the kth vector.
    real_t P[6*3];
    eval_basis_matrix(dpdx[k], dpdy[k], dpdz[k], P);

    // Compute N * D * P (a 3x3 matrix) for the kth vector.
    real_t DP[6*3];
    real_t alpha = 1.0, beta = 0.0;
    char no_trans = 'N';
    int six = 6, three = 3;
    rgemm(&no_trans, &no_trans, &six, &three, &six, &alpha, D, &six, P, &six, 
          &beta, DP, &six);
    real_t NDP[6*3];
    rgemm(&no_trans, &no_trans, &three, &three, &six, &alpha, N, &three, DP, &six,
          &beta, NDP, &three);
    I[0][0][k] = NDP[0]; I[0][1][k] = NDP[3]; I[0][2][k] = NDP[6];
    I[1][0][k] = NDP[1]; I[1][1][k] = NDP[4]; I[1][2][k] = NDP[7];
    I[2][0][k] = NDP[2]; I[2][1][k] = NDP[5]; I[2][2][k] = NDP[8];
  }
}

gmls_functional_t* elastic_gmls_functional_new(real_t E, real_t nu,
                                               int degree,
                                               point_cloud_t* points,
                                               real_t* subdomain_extents,
                                               real_t delta)
{
  ASSERT(E > 0.0);
  ASSERT(nu > 0.0);
  ASSERT(degree >= 2);
  ASSERT(degree <= 4);
  elastic_t* elastic = polymec_malloc(sizeof(elastic_t));
  elastic->E = E;
  elastic->nu = nu;
  elastic->degree = degree;
  surface_integral_t* Q = mlpg_cube_surface_integral_new(points, subdomain_extents, degree, delta);
  gmls_functional_vtable vtable = {.eval_integrands = elastic_eval_integrands,
                                   .dtor = polymec_free};
  return surface_gmls_functional_new("Elasticity Equations", elastic, vtable, 3, Q);
}

