// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_POLYNOMIAL_FIT_H
#define POLYMEC_POLYNOMIAL_FIT_H

#include "core/polymec.h"
#include "core/point.h"

// An object of this type represents a multi-component polynomial fit, 
// achieved by solving an overdetermined system of equations using least 
// squares methods.
typedef struct polynomial_fit_t polynomial_fit_t;

// This enumerated type designates methods used to solve the least-squares
// systems that arise from polynomial fits.
typedef enum
{
  QR_FACTORIZATION,
  ORTHOGONAL_FACTORIZATION,
  SINGULAR_VALUE_DECOMPOSITION,
  NORMAL_EQUATIONS
} polynomial_fit_solver_t;

// Creates a new empty least squares system for fitting multi-component 
// scatter data to a set of degree p polynomials. The type of underlying 
// solver must also be given.
polynomial_fit_t* polynomial_fit_new(int num_components, 
                                     int p,
                                     polynomial_fit_solver_t solver_type);

// Frees the polynomial fit.
void polynomial_fit_free(polynomial_fit_t* fit);

// Returns the degree of the polynomial fit.
int polynomial_fit_degree(polynomial_fit_t* fit);

// Returns the number of components in the polynomial fit.
int polynomial_fit_num_components(polynomial_fit_t* fit);

// Adds an equation to the least squares fit that attempts to interpolate 
// the scatter datum u for the given component at the point x. The equation
// is weighted with the given (positive) weight.
void polynomial_fit_add_scatter_datum(polynomial_fit_t* fit, int component, 
                                      real_t u, point_t* x, real_t weight);

// Adds an equation to the least squares system that satisfies the 
// relationship alpha * u + beta * n o grad u = gamma for the given component 
// at the point x, where alpha, beta, and gamma are real-valued quantities, 
// n is a vector, and u is the quantity being fitted to the polynomial. This 
// constraint is often referred to as a Robin boundary condition. The equation
// is weighted with the given (positive) weight.
void polynomial_fit_add_robin_bc(polynomial_fit_t* fit, int component, 
                                 real_t alpha, real_t beta, vector_t* n, real_t gamma, 
                                 point_t* x, real_t weight);

// Resets the least squares system governing the polynomial fit, removing all 
// its equations and recentering it at the point x0. If x0 is NULL, it will 
// be set to the origin.
void polynomial_fit_reset(polynomial_fit_t* fit, point_t* x0);

// Returns the number of equations in the underlying least squares system.
int polynomial_fit_num_equations(polynomial_fit_t* fit);

// Computes the polynomial fit, solving the least squares system using 
// singular value decomposition (SVD). After this is called, the fit can be 
// used to interpolate the scatter data and its derivatives.
void polynomial_fit_compute(polynomial_fit_t* fit);

// Evaluates the computed polynomial fit at the given point in space.
void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value);

// Evaluates the given partial derivative of the computed polynomial fit at 
// the given point in space.
void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               point_t* x, 
                               real_t* deriv);

// Returns the dimension of the polynomial basis for this fit. This is the 
// number of coefficients in the polynomial that represents fitted data.
int polynomial_fit_dimension(polynomial_fit_t* fit);

// Retrieves the coefficients of the polynomial computed by this object for 
// the given component, using polynomial_fit_compute(). These coefficients 
// will be stored in the given coeffs array, which should be large enough to 
// store all of the coefficients.
void polynomial_fit_get_coeffs(polynomial_fit_t* fit, 
                               int component,
                               real_t* coeffs);

// Writes a text representation of the polynomial fit to the given stream.
void polynomial_fit_fprintf(polynomial_fit_t* fit, FILE* stream);

#endif
