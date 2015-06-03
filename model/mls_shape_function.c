// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/mls_shape_function.h"

typedef struct
{
  int i;
  int poly_degree, basis_dim;
  polynomial_t* P;
  kernel_function_t* W;
  point_cloud_t* domain;
  stencil_t* neighborhoods;
  real_t* smoothing_lengths;
  real_t* basis;
} mls_t;

static int mls_neighborhood_size(void* context, int i)
{
  mls_t* mls = context;
  return stencil_size(mls->neighborhoods, i);
}

static void mls_get_neighborhood_points(void* context, int i, point_t* points)
{
  mls_t* mls = context;
  int pos = 0, j, k = 0;
  points[k++] = mls->domain->points[i];
  while (stencil_next(mls->neighborhoods, i, &pos, &j, NULL))
    points[k++] = mls->domain->points[j];
}

static void mls_set_neighborhood(void* context, int i)
{
  mls_t* mls = context;

  // Compute the basis vectors.
  int num_points = stencil_size(mls->neighborhoods, i);
  mls->basis = polymec_realloc(mls->basis, sizeof(real_t) * mls->basis_dim * num_points);
  int pos = 0, j;
  while (stencil_next(mls->neighborhoods, i, &pos, &j, NULL))
    polynomial_compute_basis(mls->poly_degree, 0, 0, 0, &mls->domain->points[j], &mls->basis[mls->basis_dim*j]);
}

static void mls_compute(void* context, 
                        int i, 
                        point_t* x,
                        real_t* values, 
                        vector_t* gradients)
{
}

static void mls_dtor(void* context)
{
  mls_t* mls = context;
  mls->P = NULL;
  polymec_free(mls->basis);
  polymec_free(mls);
}

shape_function_t* mls_shape_function_new(int polynomial_degree,
                                         kernel_function_t* W,
                                         point_cloud_t* domain,
                                         stencil_t* neighborhoods,
                                         real_t* smoothing_lengths)
{
  ASSERT(polynomial_degree >= 0);
  ASSERT(polynomial_degree <= 4);

  mls_t* mls = polymec_malloc(sizeof(mls_t));
  mls->W = W;

  // Set up a polynomial to evaluate the moment matrix.
  mls->poly_degree = polynomial_degree;
  mls->basis_dim = polynomial_basis_dim(polynomial_degree);
  real_t poly_coeffs[mls->basis_dim];
  for (int i = 0; i < mls->basis_dim; ++i)
    poly_coeffs[i] = 1.0;
  mls->P = polynomial_new(polynomial_degree, poly_coeffs, NULL);
  mls->domain = domain;
  mls->neighborhoods = neighborhoods;
  mls->smoothing_lengths = smoothing_lengths;
  mls->basis = NULL;
  shape_function_vtable vtable = {.neighborhood_size = mls_neighborhood_size,
                                  .get_neighborhood_points = mls_get_neighborhood_points,
                                  .compute = mls_compute,
                                  .dtor = mls_dtor};
  char name[1024];
  snprintf(name, 1023, "MLS shape function (p = %d)", polynomial_degree);
  return shape_function_new(name, mls, vtable);
}

