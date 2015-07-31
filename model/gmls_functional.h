// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GMLS_FUNCTIONAL_H
#define POLYMEC_GMLS_FUNCTIONAL_H

#include "model/point_weight_function.h"

// This class represents a set of functionals {lambda(u)} that represent weak 
// forms involving spatial integrals of some expression of a solution u and or 
// its derivatives, projected to a polynomial vector space. The approximation 
// of the functionals by linear combinations of basis vectors is an extension 
// of the Moving Least Squares (MLS) approximation that is referred to as 
// Generalized Moving Least Squares (see Mirzaei et al, "On generalized moving
// least squares and diffuse derivatives", IMA J. Numer. Anal. 32, 983 (2012)).
// This is useful in constructing Meshless Local Petrov-Galerkin (MLPG) methods
// that do not rely on shape functions, such as the Direct MLPG (DMLPG) method.
typedef struct gmls_functional_t gmls_functional_t;

// This virtual table defines the behavior of a GMLS functional.
typedef struct
{
  // This function returns the number of nodes contributing to the ith 
  // subdomain.
  int (*num_nodes)(void* context, int i);

  // This function retrieves the indices of the nodes contributing to the 
  // ith subdomain, placing them into the nodes array.
  void (*get_nodes)(void* context, int i, int* nodes);

  // This function retrieves the locations of the nodes with the given indices.
  void (*get_points)(void* context, int* nodes, int num_nodes, point_t* points);

  // This function returns the number of quadrature points for the quadrature 
  // rule underlying this functional on the ith subdomain.
  int (*num_quad_points)(void* context, int i);

  // This function computes the N quadrature points and weights for the ith 
  // subdomain on which the functional is to be evaluated.
  void (*get_quadrature)(void* context, int i, int N, point_t* points, real_t* weights);

  // This function computes the polynomial basis vectors composing the 
  // terms in the functional, evaluating them at the point x. These basis 
  // vectors are linear combinations of the standard polynomial basis vectors 
  // for polynomials of the given degree and their derivatives, all of which 
  // can be computed using the polynomial_compute_basis() function. The ith 
  // basis vector is stored in P[i].
  void (*compute_basis)(void* context, int degree, point_t* x, real_t* P);

  // This function evaluates the integrand of the functional applied to each 
  // of the Q polynomial basis vectors passed into the P array (where P[i] is 
  // the only nonzero term in the ith basis vector), evaluated at the given 
  // point x at time t. The values of the integrands are stored in the 
  // integrands array.
  void (*eval_integrands)(void* context, real_t t, point_t* x, real_t* P, real_t* integrands);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); // Destructor
} gmls_functional_vtable;

// Creates a generalized MLS functional with the given name, polynomial degree, 
// weight function, stencil, context, and virtual table. The weight function 
// is consumed by the functional.
gmls_functional_t* gmls_functional_new(const char* name,
                                       int poly_degree,
                                       point_weight_function_t* W,
                                       void* context,
                                       gmls_functional_vtable vtable);

// Destroys the given GMLS functional.
void gmls_functional_free(gmls_functional_t* functional);

// Returns the number of nodes contributing to the functional on the ith 
// subdomain.
int gmls_functional_num_nodes(gmls_functional_t* functional, int i);

// Evaluates the GMLS coefficients that directly approximate the functional 
// on the ith subdomain at time t. The values of these coefficients are placed 
// in the coeffs array, which should be equal in size to the number of nodes 
// in the ith subdomain. The indices of the nodes contributing to the 
// subdomain are stored in the nodes array.
void gmls_functional_compute_coeffs(gmls_functional_t* functional,
                                    real_t t,
                                    int i,
                                    real_t* coeffs, 
                                    int* nodes);

#endif
