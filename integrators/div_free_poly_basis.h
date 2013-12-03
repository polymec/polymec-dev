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
// Yes, it is that specific and esoteric, and its name is abbreviated 
// to something confusing. If you need this, you'll be glad it exists. 
// Otherwise, just ignore it. :-) Objects of this type are garbage-collected.
typedef struct div_free_poly_basis_t div_free_poly_basis_t;

// Constructs a divergence-free polynomial basis suitable for representing 
// functions that can be exactly represented by polynomials of the given 
// degree. This basis is obtained by performing a Gram-Schmidt orthogonalization
// on the naive basis whose vector components are monomials.
div_free_poly_basis_t* div_free_poly_basis_new(int degree);

// Returns the dimension of the given divergence-free polynomial basis.
int div_free_poly_basis_dim(div_free_poly_basis_t* basis);

// Allows iteration over the polynomials in the basis. Returns true 
// if terms remain, false if this iteration yields nothing. Gives access to 
// polynomials representing x, y, and z components of each vector within the 
// basis. Set pos to 0 to reset iteration.
bool div_free_poly_basis_next(div_free_poly_basis_t* basis, 
                              int* pos, 
                              polynomial_t** x_polynomial,
                              polynomial_t** y_polynomial,
                              polynomial_t** z_polynomial);

#endif
