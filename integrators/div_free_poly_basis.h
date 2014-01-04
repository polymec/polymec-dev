// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
