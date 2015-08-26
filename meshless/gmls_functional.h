// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GMLS_FUNCTIONAL_H
#define POLYMEC_GMLS_FUNCTIONAL_H

#include "core/polynomial.h"
#include "integrators/volume_integral.h"
#include "integrators/surface_integral.h"

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
  // This function evaluates the integrand of the functional applied to each 
  // of the numcomp * Q polynomial basis vectors represented by the given 
  // gmls_poly_basis object, to be evaluated at the given point x at time t. 
  // If the weak form is integrated over the boundary of the subdomain, n will 
  // point to the normal vector corresponding to x; otherwise it will be NULL. 
  // The values of the solution may be given in the component-minor-ordered 
  // array solution; if not given, solution will be NULL. The values of the 
  // integrands are stored in the integrands array, which should be sized to 
  // Nc * Nc * dim, where Nc is the number of solution components and dim 
  // is the (maximum) dimension of the polynomial basis, since the weak form 
  // for the given component can involve terms from all solution components. 
  // The integrands are stored in component-minor order, meaning that they 
  // can be interpreted as a 3D array I in which I[i][j][k] refers to the 
  // coefficient of the kth basis vector for the jth component of the weak 
  // form corresponding to the the ith component of the solution.
  void (*eval_integrands)(void* context, real_t t, 
                          multicomp_poly_basis_t* basis, 
                          point_t* x, vector_t* n, real_t* solution,
                          real_t* integrands);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); // Destructor
} gmls_functional_vtable;

// Creates a generalized MLS functional with the given name, context, virtual 
// table, number of components, and volume integral rule.
// The solution is assumed to be a vector of degrees of freedom 
// expressed in node-major order.
gmls_functional_t* volume_gmls_functional_new(const char* name,
                                              void* context,
                                              gmls_functional_vtable vtable,
                                              int num_components,
                                              volume_integral_t* quad_rule);

// Creates a generalized MLS functional with the given name, context, virtual 
// table, number of components, and surface integral rule.
// The solution is assumed to be a vector of degrees of freedom 
// expressed in node-major order.
gmls_functional_t* surface_gmls_functional_new(const char* name,
                                               void* context,
                                               gmls_functional_vtable vtable,
                                               int num_components,
                                               surface_integral_t* quad_rule);
 
// Destroys the given GMLS functional.
void gmls_functional_free(gmls_functional_t* functional);

// Returns the number of solution components in the functional.
int gmls_functional_num_components(gmls_functional_t* functional);

// Evaluates the functionals {lambda_j} consisting of the weak forms applied
// to the polynomials within the polynomial basis on the ith subdomain, 
// evaluated at time t. The value of the solution may be given as an array in 
// component-minor form (for nonlinear problems), or may be NULL (for linear 
// ones). The functionals are computed using the components (vectors) in the 
// given multi-component polynomial basis. The values are placed in the 
// lambdas array in component-minor form, which should be equal in length to 
// Nc * Nc * dim, where Nc is the number of solution components and dim 
// is the (maximum) dimension of the polynomial basis.
void gmls_functional_compute(gmls_functional_t* functional,
                             int i,
                             real_t t,
                             multicomp_poly_basis_t* poly_basis,
                             real_t* solution,
                             real_t* lambdas);

// Evaluates the integrands applied to the polynomials within the polynomial 
// basis at time t on the point x, with the normal vector n (if the functional 
// is defined at the boundary of the subdomain). The value of the solution may 
// be given as an array in component-minor form (for nonlinear problems), or 
// may be NULL (for linear ones). The integrands are placed in the integrands 
// array, which should be equal in length to Nc * Nc * dim, where Nc is the 
// number of solution components and dim is the (maximum) dimension of the 
// polynomial basis. The integrands are stored in component-minor order, 
// meaning that they can be interpreted as a 3D array I in which I[i][j][k] 
// refers to the coefficient of the kth basis vector for the jth component of 
// the weak form corresponding to the the ith component of the solution.
void gmls_functional_eval_integrands(gmls_functional_t* functional,
                                     real_t t,
                                     multicomp_poly_basis_t* poly_basis,
                                     point_t* x, vector_t* n, 
                                     real_t* solution,
                                     real_t* integrands);

#endif
