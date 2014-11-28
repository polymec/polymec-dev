// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#ifndef POLYMEC_polynomial_fit_H
#define POLYMEC_polynomial_fit_H

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
  SINGULAR_VALUE_DECOMPOSITION
} polynomial_fit_solver_t;

// Creates a new empty least squares system for fitting multi-component 
// scatter data to a set of degree p polynomials. The type of underlying 
// solver must also be given.
polynomial_fit_t* polynomial_fit_new(int num_components, 
                                     int p,
                                     polynomial_fit_solver_t solver_type);

// Frees the polynomial fit.
void polynomial_fit_free(polynomial_fit_t* fit);

// Sets the polynomial fit to use uniform weighting. This is the default
// setting for polynomial fits.
void polynomial_fit_set_unweighted(polynomial_fit_t* fit);

// Sets up an order-p spline weighting function to use for the next 
// polynomial fit, in the form 
//
// W(x, x0) = { sum(i, ai * (|x - x0|/h)**i), |x - x0| <= h,
//            { 0                           , |x - x0| >  h
// 
// where i runs from 0 to p, h is the radius of compact support, and 
// {ai} are coefficients satisfying boundary conditions that cause W and its 
// derivatives to vanish at the spherical boundary |x - x0| == h. Orders 1 
// through 4 are supported.
void polynomial_fit_set_spline_weights(polynomial_fit_t* fit, 
                                       int p,
                                       real_t h);

// Sets up an inverse distance weighting function to use for the next 
// polynomial fit, in the form 
//                        W0
// W(x, x0) = ----------------------------,
//            ((|x - x0|/h)**p + epsilon**p)
// 
// where W0 is a constant, epsilon is a "softening parameter," h is a 
// scaling parameter, and p is the exponent dictating the power law.
void polynomial_fit_set_inverse_power_weights(polynomial_fit_t* fit, 
                                              real_t W0, 
                                              real_t p,
                                              real_t h,
                                              real_t epsilon);

// Returns the degree of the polynomial fit.
int polynomial_fit_degree(polynomial_fit_t* fit);

// Returns the number of components in the polynomial fit.
int polynomial_fit_num_components(polynomial_fit_t* fit);

// Adds an equation to the least squares fit that attempts to interpolate 
// the scatter datum u for the given component at the point x.
void polynomial_fit_add_scatter_datum(polynomial_fit_t* fit, int component, 
                                      real_t u, point_t* x);

// Adds an equation to the least squares system that satisfies the 
// relationship alpha * u + beta * n o grad u = gamma for the given component 
// at the point x, where alpha, beta, and gamma are real-valued quantities, 
// n is a vector, and u is the quantity being fitted to the polynomial. This 
// constraint is often referred to as a Robin boundary condition.
void polynomial_fit_add_robin_bc(polynomial_fit_t* fit, int component, 
                                 real_t alpha, real_t beta, vector_t* n, real_t gamma, 
                                 point_t* x);

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

// Writes a text representation of the polynomial fit to the given stream.
void polynomial_fit_fprintf(polynomial_fit_t* fit, FILE* stream);

#endif
