// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GMLS_FUNCTIONAL_H
#define POLYMEC_GMLS_FUNCTIONAL_H

#include "core/polynomial.h"

// This class represents a functional lambda(u) that represents a weak 
// form for a solution u and or its derivatives, projected to a polynomial 
// vector space. The approximation of the functional by linear combinations 
// of solution data on nodes is an extension of the Moving Least Squares (MLS) 
// approximation that is referred to as Generalized Moving Least Squares 
// (see Mirzaei et al, "On generalized moving least squares and diffuse 
// derivatives", IMA J. Numer. Anal. 32, 983 (2012)).
//
// This construct is useful in constructing Meshless Local Petrov-Galerkin 
// (MLPG) methods that do not rely on shape functions, such as the Direct 
// MLPG (DMLPG) method.
typedef struct gmls_functional_t gmls_functional_t;

// This virtual table defines the behavior of a GMLS functional.
typedef struct
{
  // This function returns the number of quadrature points for the quadrature 
  // rule underlying this functional on the ith subdomain.
  int (*num_quad_points)(void* context, int i);

  // This function computes the N quadrature points and weights for the ith 
  // subdomain on which the functional is to be evaluated.
  void (*get_quadrature)(void* context, int i, int N, point_t* points, real_t* weights);

  // This function evaluates the integrand of the functional applied to each 
  // of the numcomp * Q polynomial basis vectors represented by the given 
  // gmls_poly_basis object, to be evaluated for the given component at the 
  // given point x at time t. The values of the integrands are stored in the 
  // integrands array, which should be sized to the num_comps * dimension of 
  // the given component of the polynomial basis, since the weak form for the 
  // given component can involve terms from all solution components. The 
  // integrands are stored in basis-major (component-minor) order.
  void (*eval_integrands)(void* context, int component, real_t t, point_t* x, 
                          multicomp_poly_basis_t* basis, real_t* integrands);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); // Destructor
} gmls_functional_vtable;

// Creates a generalized MLS functional with the given name, context, virtual 
// table, (multicomponent) polynomial basis, and weight function. The 
// weight function is consumed by the functional. The solution is assumed to 
// be a vector of degrees of freedom expressed in node-major order.
gmls_functional_t* gmls_functional_new(const char* name,
                                       void* context,
                                       gmls_functional_vtable vtable,
                                       multicomp_poly_basis_t* poly_basis);

// Destroys the given GMLS functional.
void gmls_functional_free(gmls_functional_t* functional);

// Returns the multi-component polynomial basis associated with this functional.
multicomp_poly_basis_t* gmls_functional_basis(gmls_functional_t* functional);

// Returns the number of solution components in the functional.
int gmls_functional_num_components(gmls_functional_t* functional);

// Evaluates the functionals {lambda_j} consisting of the weak forms applied
// to the polynomials within the polynomial basis for the given component
// on the ith subdomain evaluated at time t. The functionals are placed in 
// the lambdas array, which should be equal in length to the number of 
// components * the dimension of the polynomial basis, filled in basis-major 
// (component-minor) order.
void gmls_functional_compute(gmls_functional_t* functional,
                             int component,
                             int i,
                             real_t t,
                             real_t* lambdas);

#endif
