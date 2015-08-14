// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "meshless/mlpg_quadrature.h"
#include "poisson_gmls_functional.h"

typedef struct
{
  int degree;
  point_cloud_t* points;
  real_t* subdomain_extents;
} poisson_t;

static int poisson_num_quad_points(void* context, int i)
{
  poisson_t* poisson = context;
  return 6 * pow(poisson->degree-1, 3);
}

static void poisson_get_quadrature_p2(void* context, int i, int N, 
                                      point_t* points, real_t* weights, vector_t* normals)
{
  poisson_t* poisson = context;
  point_t xi = poisson->points->points[i];
  real_t hi = poisson->subdomain_extents[i];
  memset(normals, 0, sizeof(vector_t) * N);
  points[0] = xi;
  points[0].x = xi.x - 0.5*hi;
  weights[0] = 1.0;
  normals[0].x = -1.0;
  points[1] = xi;
  points[1].x = xi.x + 0.5*hi;
  weights[1] = 1.0;
  normals[1].x = 1.0;
  points[2] = xi;
  points[2].y = xi.y - 0.5*hi;
  weights[2] = 1.0;
  normals[2].y = -1.0;
  points[3] = xi;
  points[3].y = xi.y + 0.5*hi;
  weights[3] = 1.0;
  normals[3].y = 1.0;
  points[4] = xi;
  points[4].z = xi.z - 0.5*hi;
  weights[4] = 1.0;
  normals[4].y = -1.0;
  points[5] = xi;
  points[5].z = xi.z + 0.5*hi;
  weights[5] = 1.0;
  normals[5].y = 1.0;
}

static void poisson_eval_integrands(void* context, int component, 
                                    real_t t, point_t* x, vector_t* n,
                                    multicomp_poly_basis_t* basis, 
                                    real_t* integrands)
{
  int dim = multicomp_poly_basis_dim(basis);
  real_t dpdx[dim], dpdy[dim], dpdz[dim];
  multicomp_poly_basis_compute(basis, 0, 1, 0, 0, x, dpdx);
  multicomp_poly_basis_compute(basis, 0, 0, 1, 0, x, dpdy);
  multicomp_poly_basis_compute(basis, 0, 0, 0, 1, x, dpdz);
  for (int i = 0; i < dim; ++i)
    integrands[i] = dpdx[i] * n->x + dpdy[i] * n->y + dpdz[i] * n->z;
}

gmls_functional_t* poisson_gmls_functional_new(int degree,
                                               point_cloud_t* points,
                                               real_t* subdomain_extents,
                                               real_t delta)
{
  ASSERT(degree >= 2);
  ASSERT(degree <= 4);
  poisson_t* poisson = polymec_malloc(sizeof(poisson_t));
  poisson->degree = degree;
  multicomp_poly_basis_t* P = standard_multicomp_poly_basis_new(1, degree);
  surface_integral_t* Q = mlpg_cube_surface_integral_new(points, subdomain_extents, degree, delta);
  gmls_functional_vtable vtable = {.eval_integrands = poisson_eval_integrands,
                                   .dtor = polymec_free};
  return surface_gmls_functional_new("Poisson's Equation", poisson, vtable, P, Q);
}

