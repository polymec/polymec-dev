// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

// This fills the given array with the components of the given partial 
// derivative (of the standard polynomial basis of the given degree, 
// evaluated at the point x. 
void polynomial_compute_basis(int degree, 
                              int x_deriv, int y_deriv, int z_deriv, 
                              point_t* x, 
                              real_t* basis);

// Creates a polynomial of the given degree from the given set of coefficients,
// expanded about the given point x0. If x0 is NULL, x0 = 0.
// The coefficients correspond to monomial terms in the standard 3D polynomial 
// basis, moving from across the corresponding row in Pascal's 
// (hyper-)triangle. coeffs is an array of size polynomial_basis_dim(degree), 
// and data is copied from coeffs into the polynomial object.
polynomial_t* polynomial_new(int degree, real_t* coeffs, point_t* x0);

// Creates a polynomial of the given degree using the given monomials
// expressed in terms of coefficients and powers of x, y, and z. The 
// polynomial is expanded about x0, unless x0 is NULL, in which case it is 
// expanded about the origin.
polynomial_t* polynomial_from_monomials(int degree, int num_monomials, real_t* coeffs, 
                                        int* x_powers, int* y_powers, int* z_powers, 
                                        point_t* x0);

// Clones an existing polynomial p.
polynomial_t* polynomial_clone(polynomial_t* p);

// Creates a copy of the given polynomial scaled by the given factor.
polynomial_t* scaled_polynomial_new(polynomial_t* p, real_t factor);

// Returns the degree of the polynomial.
int polynomial_degree(polynomial_t* p);

// Returns the number of terms in the polynomial.
int polynomial_num_terms(polynomial_t* p);

// Returns an internal array of the polynomial's coefficients. Can be used 
// to get or set the coefficients.
real_t* polynomial_coeffs(polynomial_t* p);

// Returns the point about which the polynomial is expanded.
point_t* polynomial_x0(polynomial_t* p);

// Returns the scale factor (if any) associated with the polynomial.
real_t polynomial_scale_factor(polynomial_t* p);

// "Shifts" the polynomial so that its origin is x0.
void polynomial_shift(polynomial_t* p, point_t* x0);

// Scales the polynomial by the given scale factor.
void polynomial_scale(polynomial_t* p, real_t factor);

// Evaluates the polynomial at the given point x.
real_t polynomial_value(polynomial_t* p, point_t* x);

// Evaluates the (mixed partial) derivative of the polynomial at the given point x.
real_t polynomial_deriv_value(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x);

// Allows iteration over the (monomial) terms of the polynomial. Returns true 
// if terms remain, false if this iteration yields nothing. Set pos to 0 to 
// reset iteration.
bool polynomial_next(polynomial_t* p, int* pos, real_t* coeff, int* x_power, int* y_power, int* z_power);

// Adds the polynomial q (times the given factor) to the polynomial p in-place.
void polynomial_add(polynomial_t* p, real_t factor, polynomial_t* q);

// Returns a newly-created polynomial that is the product of the two 
// polynomials p and q.
polynomial_t* polynomial_product(polynomial_t* p, polynomial_t* q);

// Returns a newly-created polynomial that is the partial derivative of 
// the given polynomial with the given degrees in x, y, and z.
polynomial_t* polynomial_derivative(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv);

// Returns the index of the coefficient within the polynomial whose powers 
// match the given x_power, y_power, and z_power. If such a term is not 
// found in the polynomial, returns -1.
int polynomial_index(polynomial_t* p, int x_power, int y_power, int z_power);

// Returns true if the two given polynomials are equal, false if not.
bool polynomial_equals(polynomial_t* p, polynomial_t* q);

// Returns an sp_func corresponding to the given polynomial.
sp_func_t* polynomial_sp_func(polynomial_t* p);

// Writes a text representation of the polynomial to the given stream.
void polynomial_fprintf(polynomial_t* p, FILE* stream);

// This is a simple class that represents a non-standard polynomial basis.
// Objects of this type are garbage-collected.
typedef struct poly_basis_t poly_basis_t;

// Creates a new polynomial basis from the given set of polynomials.
poly_basis_t* poly_basis_new(int dimension, polynomial_t** polynomials);

// Creates a new standard polynomial basis of the given dimension.
poly_basis_t* standard_poly_basis_new(int degree);

// Returns the dimension of the polynomial basis.
int poly_basis_dim(poly_basis_t* basis);

// Returns the degree of the polynomial basis.
int poly_basis_degree(poly_basis_t* basis);

// Steps through the vectors for the polynomial basis, placing each polynomial
// into *p. Set *pos to 0 to reset the traversal. Returns true if the 
// traversal yields another vector, false if it is finished.
bool poly_basis_next(poly_basis_t* basis,
                     int* pos,
                     polynomial_t** p);

// Computes the polynomial basis, evaluated at the point x. The basis 
// values are placed in the values array (sized to match the dimension of the 
// basis). x_deriv, y_deriv, and z_deriv are the indices of the desired x, y, 
// and z derivatives for the computed basis.
void poly_basis_compute(poly_basis_t* basis, 
                        int x_deriv, int y_deriv, int z_deriv,
                        point_t* x,
                        real_t* values);

// "Shifts" the polynomial basis so that its origin is x0.
void poly_basis_shift(poly_basis_t* basis, point_t* x0);

// Scales the polynomial basis by the given scale factor.
void poly_basis_scale(poly_basis_t* basis, real_t factor);

// Returns true if the two given polynomial bases are equal, false if not.
bool poly_basis_equals(poly_basis_t* basis1, poly_basis_t* basis2);

// This is a simple class that represents a multi-component polynomial basis.
// Objects of this type are garbage-collected.
typedef struct multicomp_poly_basis_t multicomp_poly_basis_t;

// Creates a new multi-component polynomial basis from the given array of
// single-component polynomial bases.
multicomp_poly_basis_t* multicomp_poly_basis_new(int num_components, 
                                                 poly_basis_t** component_bases);

// Creates a new multi-component polynomial from a set of standard 
// component-wise polynomial bases, all of the given degree.
multicomp_poly_basis_t* standard_multicomp_poly_basis_new(int num_components,
                                                          int degree);

// Creates a new multi-component polynomial from a set of standard 
// component-wise polynomial bases, with the ith component having degree
// degrees[i].
multicomp_poly_basis_t* var_standard_multicomp_poly_basis_new(int num_components,
                                                              int* degrees);

// Returns the number of solution components in the GMLS polynomial basis.
int multicomp_poly_basis_num_comp(multicomp_poly_basis_t* basis);

// Returns the dimension of the multi-component polynomial basis, which is 
// the maximum dimension found for the underlying polynomial bases.
int multicomp_poly_basis_dim(multicomp_poly_basis_t* basis);

// Returns the degree of the multi-component polynomial basis.
int multicomp_poly_basis_degree(multicomp_poly_basis_t* basis);

// Steps through the vectors for the multi-component polynomial basis. p[i] 
// is set to a polynomial representing the present vector for component i.
// Returns true if the traversal yields another vector, false if it is 
// finished.
bool multicomp_poly_basis_next(multicomp_poly_basis_t* basis,
                               int* pos,
                               polynomial_t*** p);

// Computes the polynomial basis for the given component, evaluated at 
// the point x, storing the values in the values array (sized to match 
// the maximum dimension of the basis for all components). x_deriv, y_deriv, 
// and z_deriv are the indices of the desired x, y, and z derivatives for the 
// computed basis.
void multicomp_poly_basis_compute(multicomp_poly_basis_t* basis,
                                  int component, 
                                  int x_deriv, int y_deriv, int z_deriv,
                                  point_t* x,
                                  real_t* values);

// "Shifts" the polynomial basis so that its origin is x0.
void multicomp_poly_basis_shift(multicomp_poly_basis_t* basis, point_t* x0);

// Scales the polynomial basis by the given scale factor.
void multicomp_poly_basis_scale(multicomp_poly_basis_t* basis, real_t factor);

// Returns true if the two multi-component polynomial bases are equal, 
// false if not.
bool multicomp_poly_basis_equals(multicomp_poly_basis_t* basis1, 
                                 multicomp_poly_basis_t* basis2);

// Returns true if all of the components have equal polynomial bases, 
// false if not.
bool multicomp_poly_basis_components_equal(multicomp_poly_basis_t* basis);

#endif
