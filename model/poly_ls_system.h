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

#ifndef POLYMEC_POLY_LS_SYSTEM_H
#define POLYMEC_POLY_LS_SYSTEM_H

#include "core/polymec.h"
#include "core/point.h"

// An object of this type represents an overdetermined system of equations 
// that can be solved by fitting data to polynomials using least squares methods.
typedef struct poly_ls_system_t poly_ls_system_t;

// Creates a new empty least squares system for fitting multi-component 
// scatter data to a set of degree p polynomials.
poly_ls_system_t* poly_ls_system_new(int num_components, int p);

// Frees the least squares system.
void poly_ls_system_free(poly_ls_system_t* sys);

// Sets the center point x0 about which the underlying polynomials are centered.
// By default, these polynomials are centered at the origin.
void poly_ls_system_set_x0(poly_ls_system_t* sys, point_t* x0);

// Adds an equation to the least squares system that interpolates the datum 
// u for the given component at the point x.
void poly_ls_system_add_interpolated_datum(poly_ls_system_t* sys, int component, 
                                           real_t u, point_t* x);

// Adds an equation to the least squares system that satisfies the 
// relationship alpha * u + beta * n o grad u = gamma for the given component 
// at the point x, where alpha, beta, and gamma are real-valued quantities, 
// n is a vector, and u is the quantity being fitted to the polynomial. This 
// constraint is often referred to as a Robin boundary condition.
void poly_ls_system_add_robin_bc(poly_ls_system_t* sys, int component, 
                                 real_t alpha, real_t beta, vector_t* n, real_t gamma, 
                                 point_t* x);

// Clears the least squares system, removing all its equations.
void poly_ls_system_clear(poly_ls_system_t* sys);

// Returns the number of equations in the least squares system.
int poly_ls_system_num_equations(poly_ls_system_t* sys);

// Solves the least squares system using singular value decomposition (SVD),
// placing the solution (the coefficients of the fitting polynomials for the 
// components, stored in component-minor form) into the array x.
void poly_ls_system_solve(poly_ls_system_t* sys, real_t* x);

#endif
