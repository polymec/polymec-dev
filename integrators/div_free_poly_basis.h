// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
div_free_poly_basis_t* spherical_div_free_poly_basis_new(int degree, point_t* x0, double radius);

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
