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

  point_weight_function_t* W;
  int degree, dim;
};

gmls_functional_t* gmls_functional_new(const char* name,
                                       int poly_degree,
                                       point_weight_function_t* W,
                                       void* context,
                                       gmls_functional_vtable vtable)
{
  ASSERT(poly_degree >= 2);
  ASSERT(vtable.num_nodes != NULL);
  ASSERT(vtable.get_nodes != NULL);
  ASSERT(vtable.get_points != NULL);
  ASSERT(vtable.num_quad_points != NULL);
  ASSERT(vtable.get_quadrature != NULL);
  ASSERT(vtable.compute_basis != NULL);
  ASSERT(vtable.eval_integrands != NULL);

  gmls_functional_t* functional = polymec_malloc(sizeof(gmls_functional_t));
  functional->name = string_dup(name);
  functional->context = context;
  functional->vtable = vtable;
  functional->W = W;
  functional->degree = poly_degree;
  functional->dim = polynomial_basis_dim(poly_degree);
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

int gmls_functional_num_nodes(gmls_functional_t* functional, int i)
{
  return functional->vtable.num_nodes(functional->context, i);
}

void gmls_functional_compute_coeffs(gmls_functional_t* functional,
                                    real_t t,
                                    int i,
                                    real_t* coeffs,
                                    int* nodes)
{
  // Compute the quadrature points/weights for this subdomain.
  int num_quad_points = functional->vtable.num_quad_points(functional->context, i);
  point_t quad_points[num_quad_points];
  real_t quad_weights[num_quad_points];
  functional->vtable.get_quadrature(functional->context, i, num_quad_points, 
                             quad_points, quad_weights);

  // Loop through the points and compute the functional approximants.
  int Q = functional->dim;
  real_t lambda[Q];
  memset(lambda, 0, sizeof(real_t) * Q);
  for (int q = 0; q < num_quad_points; ++q)
  {
    point_t* xq = &quad_points[q];
    real_t wq = quad_weights[q];

    // Compute the polynomial basis for the integrand at this point.
    real_t p[Q];
    functional->vtable.compute_basis(functional->context, functional->degree, xq, p);

    // Now compute the integrands for the functional at this point.
    real_t integrands[Q];
    functional->vtable.eval_integrands(functional->context, t, xq, p, integrands);

    // Integrate.
    for (int l = 0; l < Q; ++l)
      lambda[l] += wq * integrands[l];
  }

  // Get the nodes within this subdomain.
  int num_nodes = functional->vtable.num_nodes(functional->context, i);
  functional->vtable.get_nodes(functional->context, i, nodes);
  point_t xi, xjs[num_nodes];
  functional->vtable.get_points(functional->context, &i, 1, &xi);
  functional->vtable.get_points(functional->context, nodes, num_nodes, &xi);

  // Compute the weights for the nodes.
  real_t W[num_nodes];
  for (int j = 0; j < num_nodes; ++j)
  {
    vector_t y;
    point_displacement(&xi, &xjs[j], &y);
    W[j] = point_weight_function_value(functional->W, &y);
  }

  // Compute the matrix [P]_ij = pi(xj), in column major order.
  real_t P[num_nodes*Q];
  for (int j = 0; j < num_nodes; ++j)
    functional->vtable.compute_basis(functional->context, functional->degree, &xjs[j], &P[j*Q]);

  // Compute the moment matrix A = P * W * Pt.
  real_t A[Q*Q], Ainv_lambda[Q];
  memset(A, 0, sizeof(real_t) * Q * Q);
  for (int n = 0; n < num_nodes; ++n)
  {
    for (int j = 0; j < Q; ++j)
    {
      Ainv_lambda[j] = lambda[j];
      for (int l = 0; l < Q; ++l)
        A[Q*l+j] += P[Q*n+j] * W[n] * P[num_nodes*l+n];
    }
  }

  // Factor the moment matrix.
  char uplo = 'L';
  int info;
  rpotrf(&uplo, &Q, A, &Q, &info);
  ASSERT(info == 0);

  // Compute Ainv * lambda.
  int one = 1;
  rpotrs(&uplo, &Q, &one, A, &Q, Ainv_lambda, &Q, &info);
  ASSERT(info == 0);

  // Left multiply by W * Pt.
  memset(coeffs, 0, sizeof(real_t) * num_nodes);
  for (int n = 0; n < num_nodes; ++n)
  {
    for (int j = 0; j < Q; ++j)
      coeffs[n] += P[num_nodes*j+n] * W[n] * Ainv_lambda[j];
  }
}

