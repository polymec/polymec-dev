// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GMLS_FUNCTIONAL_H
#define POLYMEC_GMLS_FUNCTIONAL_H

#include "core/polynomial.h"
#include "model/point_weight_function.h"

// This class represents a functional lambda(u) that represents a weak 
// form for a solution u and or its derivatives, projected to a polynomial 
// vector space. The approximation of the functional by linear combinations 
// of solution data on nodes is an extension of the Moving Least Squares (MLS) 
// approximation that is referred to as Generalized Moving Least Squares 
// (see Mirzaei et al, "On generalized moving least squares and diffuse 
// derivatives", IMA J. Numer. Anal. 32, 983 (2012)).
// The output of a GMLS functional is a set of coefficients {a_j} for this 
// linear combination. This construct is useful in constructing Meshless 
// Local Petrov-Galerkin (MLPG) methods that do not rely on shape functions, 
// such as the Direct MLPG (DMLPG) method.
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

  // This function evaluates the integrand of the functional applied to each 
  // of the Q polynomial basis vectors represented by the given gmls_poly_basis 
  // object, to be evaluated at the given point x at time t. The values of the 
  // integrands are stored in the integrands array.
  void (*eval_integrands)(void* context, real_t t, point_t* x, multicomp_poly_basis_t* basis, real_t* integrands);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); // Destructor
} gmls_functional_vtable;

// Creates a generalized MLS functional with the given name, context, virtual 
// table, (multicomponent) polynomial basis, and weight function. The 
// weight function is consumed by the functional.
gmls_functional_t* gmls_functional_new(const char* name,
                                       void* context,
                                       gmls_functional_vtable vtable,
                                       multicomp_poly_basis_t* poly_basis,
                                       point_weight_function_t* W);

// Destroys the given GMLS functional.
void gmls_functional_free(gmls_functional_t* functional);

// Returns the number of solution components in the functional.
int gmls_functional_num_components(gmls_functional_t* functional);

// Returns the number of nodes contributing to the functional on the ith 
// subdomain.
int gmls_functional_num_nodes(gmls_functional_t* functional, int i);

// Evaluates the GMLS coefficients that directly approximate the functional 
// on the ith subdomain at time t. The values of these coefficients are placed 
// in the coeffs array. The coeffs array is a (num_nodes * num_components) 
// array stored in node-major order, where num_nodes is the number of nodes 
// contributing to in the ith subdomain and num_components is the number of 
// solution components. The indices of the nodes contributing to the 
// subdomain are stored in the nodes array.
void gmls_functional_compute_coeffs(gmls_functional_t* functional,
                                    real_t t,
                                    int i,
                                    real_t* coeffs, 
                                    int* nodes);

#endif
