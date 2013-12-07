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

// Constructs a divergence-free basis suitable for representing 
// functions that can be exactly represented by polynomials of the given 
// degree. 
div_free_poly_basis_t* div_free_poly_basis_new(int degree);

// Returns the dimension of the given divergence-free polynomial basis.
int div_free_poly_basis_dim(div_free_poly_basis_t* basis);

// Computes an array of vectors belonging to the divergence-free polynomial 
// basis, evaluated at the given point x. vectors is an array of vectors 
// equal in length to the dimension of the basis that will store the computed 
// vectors in the basis.
void div_free_poly_basis_compute(div_free_poly_basis_t* basis,
                                 point_t* x, 
                                 vector_t** vectors);

#endif
