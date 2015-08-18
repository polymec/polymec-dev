// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "meshless/gmls_functional.h"

struct gmls_functional_t 
{
  char *name;
  void* context;
  gmls_functional_vtable vtable;

  // Quadrature rules.
  volume_integral_t* volume_quad_rule;
  surface_integral_t* surface_quad_rule;

  int num_comp;
};

gmls_functional_t* volume_gmls_functional_new(const char* name,
                                              void* context,
                                              gmls_functional_vtable vtable,
                                              int num_components,
                                              volume_integral_t* quad_rule)
            
{
  ASSERT(vtable.eval_integrands != NULL);
  ASSERT(num_components > 0);

  gmls_functional_t* functional = polymec_malloc(sizeof(gmls_functional_t));
  functional->name = string_dup(name);
  functional->context = context;
  functional->vtable = vtable;
  functional->volume_quad_rule = quad_rule;
  functional->surface_quad_rule = NULL;
  functional->num_comp = num_components;
  return functional;
}

gmls_functional_t* surface_gmls_functional_new(const char* name,
                                               void* context,
                                               gmls_functional_vtable vtable,
                                               int num_components,
                                               surface_integral_t* quad_rule)
            
{
  ASSERT(vtable.eval_integrands != NULL);
  ASSERT(num_components > 0);

  gmls_functional_t* functional = polymec_malloc(sizeof(gmls_functional_t));
  functional->name = string_dup(name);
  functional->context = context;
  functional->vtable = vtable;
  functional->num_comp = num_components;
  functional->volume_quad_rule = NULL;
  functional->surface_quad_rule = quad_rule;
  return functional;
}

void gmls_functional_free(gmls_functional_t* functional)
{
  if ((functional->context != NULL) && (functional->vtable.dtor != NULL))
    functional->vtable.dtor(functional->context);
  polymec_free(functional->name);
  polymec_free(functional);
}

int gmls_functional_num_components(gmls_functional_t* functional)
{
  return functional->num_comp;
}

static void compute_integral(gmls_functional_t* functional,
                             int component,
                             real_t t,
                             multicomp_poly_basis_t* poly_basis,
                             point_t* quad_points,
                             real_t* quad_weights,
                             vector_t* quad_normals,
                             int num_quad_points,
                             real_t* lambdas)
{
  bool on_boundary = (quad_normals != NULL);

  // Loop through the points and compute the lambda matrix of functional 
  // approximants.
  int num_comp = functional->num_comp;
  int basis_dim = multicomp_poly_basis_dim(poly_basis);
  memset(lambdas, 0, sizeof(real_t) * basis_dim);
  for (int q = 0; q < num_quad_points; ++q)
  {
    point_t* xq = &quad_points[q];
    real_t wq = quad_weights[q];
    vector_t* nq = (on_boundary) ? &quad_normals[q] : NULL;

    // Now compute the (multi-component) integrands for the functional at 
    // this point.
    real_t integrands[num_comp*basis_dim];
    functional->vtable.eval_integrands(functional->context, component, 
                                       t, poly_basis, xq, nq, integrands);

    // Integrate.
    for (int i = 0; i < basis_dim; ++i)
{
      for (int c = 0; c < num_comp; ++c)
{
        lambdas[num_comp*i+c] += wq * integrands[num_comp*i+c];
printf("lambdas[%d] -> %g\n", num_comp*i+c, lambdas[num_comp*i+c]);
}
}
  }
}

void gmls_functional_compute(gmls_functional_t* functional,
                             int component, 
                             int i,
                             real_t t,
                             multicomp_poly_basis_t* poly_basis,
                             real_t* lambdas)
{
  if (functional->surface_quad_rule != NULL)
  {
    surface_integral_set_domain(functional->surface_quad_rule, i);
    int num_quad_points = surface_integral_num_points(functional->surface_quad_rule);
    point_t quad_points[num_quad_points];
    real_t quad_weights[num_quad_points];
    vector_t quad_normals[num_quad_points];
    surface_integral_get_quadrature(functional->surface_quad_rule, quad_points, quad_weights, quad_normals);
    compute_integral(functional, component, t, poly_basis, quad_points, quad_weights, 
                     quad_normals, num_quad_points, lambdas);
  }
  else
  {
    volume_integral_set_domain(functional->volume_quad_rule, i);
    int num_quad_points = volume_integral_num_points(functional->volume_quad_rule);
    point_t quad_points[num_quad_points];
    real_t quad_weights[num_quad_points];
    volume_integral_get_quadrature(functional->volume_quad_rule, quad_points, quad_weights);
    compute_integral(functional, component, t, poly_basis, quad_points, quad_weights, 
                     NULL, num_quad_points, lambdas);
  }
}

void gmls_functional_eval_integrands(gmls_functional_t* functional,
                                     int component,
                                     real_t t,
                                     multicomp_poly_basis_t* poly_basis,
                                     point_t* x, vector_t* n,
                                     real_t* integrands)
{
  functional->vtable.eval_integrands(functional->context, component,
                                     t, poly_basis, x, n, integrands);
}
