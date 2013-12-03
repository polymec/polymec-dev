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

#ifndef POLYMEC_POLYNOMIAL_H
#define POLYMEC_POLYNOMIAL_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/sp_func.h"

// This type represents a polynomial function (in x, y, and z) of a given degree.
// Objects of this type are garbage-collected.
typedef struct polynomial_t polynomial_t;

// This returns the number of coefficients in the standard basis for 
// polynomials of a given degree.
int polynomial_basis_dim(int degree);

// Creates a polynomial of the given degree from the given set of coefficients,
// expanded about the given point x0. If x0 is NULL, x0 = 0.
// The coefficients correspond to monomial terms in the standard 3D polynomial 
// basis, moving from across the corresponding row in Pascal's 
// (hyper-)triangle. coeffs is an array of size polynomial_basis_dim(degree), 
// and data is copied from coeffs into the polynomial object.
polynomial_t* polynomial_new(int degree, double* coeffs, point_t* x0);

// Creates a polynomial of the given degree using the given monomials
// expressed in terms of coefficients and powers of x, y, and z. The 
// polynomial is expanded about x0, unless x0 is NULL, in which case it is 
// expanded about the origin.
polynomial_t* polynomial_from_monomials(int degree, int num_monomials, double* coeffs, 
                                        int* x_powers, int* y_powers, int* z_powers, 
                                        point_t* x0);

// Returns the degree of the polynomial.
int polynomial_degree(polynomial_t* p);

// Returns the number of terms in the polynomial.
int polynomial_num_terms(polynomial_t* p);

// Returns an internal array of the polynomial's coefficients. Can be used 
// to get or set the coefficients.
double* polynomial_coeffs(polynomial_t* p);

// Returns the point about which the polynomial is expanded.
point_t* polynomial_x0(polynomial_t* p);

// Evaluates the polynomial at the given point x.
double polynomial_value(polynomial_t* p, point_t* x);

// Evaluates the (mixed partial) derivative of the polynomial at the given point x.
double polynomial_deriv(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x);

// Allows iteration over the (monomial) terms of the polynomial. Returns true 
// if terms remain, false if this iteration yields nothing. Set pos to 0 to 
// reset iteration.
bool polynomial_next(polynomial_t* p, int* pos, double* coeff, int* x_power, int* y_power, int* z_power);

// Returns an sp_func corresponding to the given polynomial.
sp_func_t* polynomial_sp_func(polynomial_t* p);

#endif
