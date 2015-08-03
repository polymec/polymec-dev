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
  int dim, num_comp;
};

gmls_functional_t* gmls_functional_new(const char* name,
                                       void* context,
                                       gmls_functional_vtable vtable,
                                       multicomp_poly_basis_t* poly_basis)
{
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
  return functional;
}

void gmls_functional_free(gmls_functional_t* functional)
{
  if ((functional->context != NULL) && (functional->vtable.dtor != NULL))
    functional->vtable.dtor(functional->context);
  polymec_free(functional->name);
  polymec_free(functional);
}

multicomp_poly_basis_t* gmls_functional_basis(gmls_functional_t* functional)
{
  return functional->basis;
}

int gmls_functional_num_components(gmls_functional_t* functional)
{
  return functional->num_comp;
}

void gmls_functional_compute(gmls_functional_t* functional,
                             int component, 
                             int i,
                             real_t t,
                             real_t* lambdas)
{
  // Compute the quadrature points/weights for this subdomain.
  int num_quad_points = functional->vtable.num_quad_points(functional->context, i);
  point_t quad_points[num_quad_points];
  real_t quad_weights[num_quad_points];
  functional->vtable.get_quadrature(functional->context, i, num_quad_points, 
                                    quad_points, quad_weights);

  // Loop through the points and compute the lambda matrix of functional 
  // approximants.
  int num_comp = functional->num_comp;
  int basis_dim = functional->dim;
  memset(lambdas, 0, sizeof(real_t) * basis_dim);
  for (int q = 0; q < num_quad_points; ++q)
  {
    point_t* xq = &quad_points[q];
    real_t wq = quad_weights[q];

    // Now compute the (multi-component) integrands for the functional at 
    // this point.
    real_t integrands[num_comp*basis_dim];
    functional->vtable.eval_integrands(functional->context, component, t, xq, functional->basis, integrands);

    // Integrate.
    for (int i = 0; i < basis_dim; ++i)
      for (int c = 0; c < num_comp; ++c)
        lambdas[num_comp*i+c] += wq * integrands[num_comp*i+c];
  }
}

