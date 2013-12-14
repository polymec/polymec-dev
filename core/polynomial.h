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

// Clones an existing polynomial p.
polynomial_t* polynomial_clone(polynomial_t* p);

// Creates a copy of the given polynomial scaled by the given factor.
polynomial_t* scaled_polynomial_new(polynomial_t* p, double factor);

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
double polynomial_deriv_value(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x);

// Allows iteration over the (monomial) terms of the polynomial. Returns true 
// if terms remain, false if this iteration yields nothing. Set pos to 0 to 
// reset iteration.
bool polynomial_next(polynomial_t* p, int* pos, double* coeff, int* x_power, int* y_power, int* z_power);

// Adds the polynomial q (times the given factor) to the polynomial p in-place.
void polynomial_add(polynomial_t* p, double factor, polynomial_t* q);

// Returns a newly-created polynomial that is the product of the two 
// polynomials p and q.
polynomial_t* polynomial_product(polynomial_t* p, polynomial_t* q);

// Returns a newly-created polynomial that is the partial derivative of 
// the given polynomial with the given degrees in x, y, and z.
polynomial_t* polynomial_derivative(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv);

// Returns an sp_func corresponding to the given polynomial.
sp_func_t* polynomial_sp_func(polynomial_t* p);

// Writes a text representation of the polynomial to the given stream.
void polynomial_fprintf(polynomial_t* p, FILE* stream);

#endif
