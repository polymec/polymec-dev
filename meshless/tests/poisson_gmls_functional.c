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

static void poisson_eval_integrands(void* context, real_t t, 
                                    multicomp_poly_basis_t* basis, 
                                    point_t* x, vector_t* n, real_t* solution,
                                    real_t* integrands)
{
  int dim = multicomp_poly_basis_dim(basis);
  real_t dpdx[dim], dpdy[dim], dpdz[dim];
  multicomp_poly_basis_compute(basis, 0, 1, 0, 0, x, dpdx);
  multicomp_poly_basis_compute(basis, 0, 0, 1, 0, x, dpdy);
  multicomp_poly_basis_compute(basis, 0, 0, 0, 1, x, dpdz);

  memset(integrands, 0, sizeof(int) * dim);
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
  surface_integral_t* Q = mlpg_cube_surface_integral_new(points, subdomain_extents, degree, delta);
  gmls_functional_vtable vtable = {.eval_integrands = poisson_eval_integrands,
                                   .dtor = polymec_free};
  return surface_gmls_functional_new("Poisson's Equation", poisson, vtable, 1, Q);
}

