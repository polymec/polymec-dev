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

#ifndef POLYMEC_LEAST_SQUARES_H
#define POLYMEC_LEAST_SQUARES_H

#include "core/polymec.h"
#include "core/point.h"

// These functions and types are intended to help with the construction of 
// least-squares fits to scattered data.

// Computes the coefficients A and b for a linear regression y = A*x + B 
// given arrays of (scalar) x and y values. Also computes the variance sigma.
void linear_regression(real_t* x, real_t* y, int N, real_t* A, real_t* B, real_t* sigma);

// This type is a weight function used for least-squares systems. Objects 
// of this type are garbage-collected.
typedef struct ls_weight_func_t ls_weight_func_t;

// This is the signature for least-squares weight functions. Arguments:
// - A context pointer for this weighting function.
// - The displacement y = x - x0 of the evaluated point from the function's "center".
// - A pointer to storage for the weighting function value.
// - A pointer to storage for the weighting function gradient.
typedef void (*ls_weight_func_eval_func)(void* context, vector_t* y, real_t* W, vector_t* grad_W);

// This is the signature for setting the "domain" (set of points) that will 
// define the extent of a least-squares weight function.
typedef void (*ls_weight_func_set_domain_func)(void* context, point_t* x0, point_t* points, int num_points);

// Destructor function for weight function context.
typedef void (*ls_weight_func_dtor)(void*);

// This virtual table must be implemented for a least-squares weight function.
typedef struct
{
  ls_weight_func_set_domain_func set_domain;
  ls_weight_func_eval_func       eval;
  ls_weight_func_dtor            dtor; 
} ls_weight_func_vtable;

// Creates a least-squares weight function from a name, a context, and a 
// virtual table.
ls_weight_func_t* ls_weight_func_new(const char* name,
                                     void* context,
                                     ls_weight_func_vtable vtable);

// Returns the name of the weight function.
const char* ls_weight_func_name(ls_weight_func_t* W);

// Sets the domain (the center point x0 and the points on which data exist) of the weight function.
void ls_weight_func_set_domain(ls_weight_func_t* W, point_t* x0, point_t* points, int num_points);

// Evaluates the value and the gradient of the weight function at the given point x.
void ls_weight_func_eval(ls_weight_func_t* W, point_t* x, real_t* value, vector_t* gradient);

// Returns the current center point x0 of the weight function. 
point_t* ls_weight_func_x0(ls_weight_func_t* W);

// An object of this type represents an equation in an overdetermined system.
// Objects of this type are garbage-collected.
typedef struct poly_ls_equation_t poly_ls_equation_t;

// Creates a new empty equation for a least squares system, appropriate for 
// fitting data to an order p polynomial centered about the point x0.
poly_ls_equation_t* poly_ls_equation_new(int p, point_t* x0);

// An object of this type represents an overdetermined system of equations 
// that can be solved by fitting data to polynomials using least squares methods.
typedef struct poly_ls_system_t poly_ls_system_t;

// Creates a new empty least squares system.
poly_ls_system_t* poly_ls_system_new();

// Frees the least squares system.
void poly_ls_system_free(poly_ls_system_t* sys);

// Adds an equation to the least squares system.
void poly_ls_system_add_equation(poly_ls_system_t* sys, poly_ls_equation_t* eq);

// Clears the least squares system, removing all its equations.
void poly_ls_system_clear(poly_ls_system_t* sys);

// Returns the number of equations in the least squares system.
int poly_ls_system_num_equations(poly_ls_system_t* sys);

// Solves the least squares system using singular value decomposition (SVD),
// placing the solution (the coefficients of the fitting polynomial) into 
// the array x.
void poly_ls_system_solve(poly_ls_system_t* sys, real_t* x);

// Computes the least squares system for a pth-order polynomial fit to a set 
// of scattered point data, centered about the point x0. This computes the 
// moment matrix and the right-hand side vector for a polynomial fit to the 
// given (scalar) data on the given points.
void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points, 
                            real_t* data, real_t* moment_matrix, real_t* rhs);

// This computes the weighted least squares system for the weighting function W, centered 
// about the point x0 and with spatial extent h.
void compute_weighted_poly_ls_system(int p, ls_weight_func_t* W, point_t* x0, 
                                     point_t* points, int num_points, real_t* data, 
                                     real_t* moment_matrix, real_t* rhs);

#endif
