// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DIV_FREE_POLY_BASIS_H
#define POLYMEC_DIV_FREE_POLY_BASIS_H

#include "core/polynomial.h"

// This type represents a basis for a vector space whose elements are 
// (vector-valued) divergence-free vectors with polynomial components.
// The polynomials are orthogonal in the sense that an inner product <u, v>
// is defined as the dot product of the polynomial vectors u and v, integrated
// over the surface of a given polytope.
typedef struct div_free_poly_basis_t div_free_poly_basis_t;

// Constructs a divergence-free basis of the given degree over the sphere
// at the given center point x0 and of the given radius.
div_free_poly_basis_t* spherical_div_free_poly_basis_new(int degree, point_t* x0, real_t radius);

// Returns the dimension of the given divergence-free polynomial basis.
int div_free_poly_basis_dim(div_free_poly_basis_t* basis);

// Steps through the vectors (with polynomial components) in the 
// divergence-free basis, in the same manner as polynomial_next().
// x, y, and z are set to polynomials representing the x, y, and z 
// components of the present vector. Returns true if the traversal yields
// another vector.
bool div_free_poly_basis_next(div_free_poly_basis_t* basis,
                              int* pos,
                              polynomial_t** x,
                              polynomial_t** y,
                              polynomial_t** z);

// Computes an array of vectors belonging to the divergence-free polynomial 
// basis, evaluated at the given point x. vectors is an array of vectors 
// equal in length to the dimension of the basis that will store the computed 
// vectors in the basis.
void div_free_poly_basis_compute(div_free_poly_basis_t* basis,
                                 point_t* x, 
                                 vector_t* vectors);

#endif
