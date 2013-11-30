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

#ifndef POLYMEC_LEAST_SQUARES_H
#define POLYMEC_LEAST_SQUARES_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/polynomial.h"

// These functions and types are intended to help with the construction of 
// least-squares fits to scattered data.

// Computes the coefficients A and b for a linear regression y = A*x + B 
// given arrays of (scalar) x and y values. Also computes the variance sigma.
void linear_regression(double* x, double* y, int N, double* A, double* B, double* sigma);

// This type is a weight function used for least-squares systems. Objects 
// of this type are garbage-collected.
typedef struct ls_weight_func_t ls_weight_func_t;

// This is the signature for least-squares weight functions. Arguments:
// - A context pointer for this weighting function.
// - The displacement y = x - x0 of the evaluated point from the function's "center".
// - A pointer to storage for the weighting function value.
// - A pointer to storage for the weighting function gradient.
typedef void (*ls_weight_func_eval_func)(void* context, vector_t* y, double* W, vector_t* grad_W);

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
void ls_weight_func_eval(ls_weight_func_t* W, point_t* x, double* value, vector_t* gradient);

// Returns the current center point x0 of the weight function. 
point_t* ls_weight_func_x0(ls_weight_func_t* W);

// Computes the least squares system for a pth-order polynomial fit to a set 
// of scattered point data, centered about the point x0. This computes the 
// moment matrix and the right-hand side vector for a polynomial fit to the 
// given (scalar) data on the given points.
void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points, 
                            double* data, double* moment_matrix, double* rhs);

// This computes the weighted least squares system for the weighting function W, centered 
// about the point x0 and with spatial extent h.
void compute_weighted_poly_ls_system(int p, ls_weight_func_t* W, point_t* x0, 
                                     point_t* points, int num_points, double* data, 
                                     double* moment_matrix, double* rhs);

// Computes the polynomial for a degree-p weighted least-squares fit with 
// the given weight function and scatter point data.
polynomial_t* ls_polynomial_new(int p, ls_weight_func_t* W, point_t* x0, 
                                point_t* points, int num_points, double* data);

// This class represents a shape function basis for a polynomial least-squares 
// fit. A shape function maps a set of data (associated a given set of points in 
// space) to a value interpolated at a given point.
// Objects of this type are garbage-collected.
typedef struct poly_ls_shape_t poly_ls_shape_t;

// Create a new shape function for a polynomial least-squares fit of order p, 
// with a weighting function, a context pointer, and its destructor.
// Set compute_gradients to true to allow the calculation of gradients of 
// shape functions at the cost of additional work in poly_ls_shape_set_domain.
// If W is NULL, an unweighted least-squares procedure will be used.
poly_ls_shape_t* poly_ls_shape_new(int p, ls_weight_func_t* W, bool compute_gradients);

// Sets the domain of the shape function: its origin x0, and its support points.
void poly_ls_shape_set_domain(poly_ls_shape_t* N, point_t* x0, point_t* points, int num_points);

// Computes the shape function basis evaluating each shape function at the point x.
// poly_ls_shape_set_domain must have been called previously.
void poly_ls_shape_compute(poly_ls_shape_t* N, point_t* x, double* values);

// Computes the gradients of the shape function basis (expanded about the 
// point x0 and fitted to the given points), evaluating each gradient at the point x.
void poly_ls_shape_compute_gradients(poly_ls_shape_t* N, point_t* x, double* values, vector_t* gradients);

// Computes the matrix and vector for the affine transformation that maps 
// the set of all solution values {phi} to a set of constrained solutions 
// {phic}. Together, the unconstrained and constrained solutions make up the 
// support of the least-squares polynomial shape function. The constrained 
// solution values obey the constraint:
//
// ai * phii + bi * dphii/dx + ci * dphii/dy + di * dphii/dz = ei
//
// at each of the points {xi} at which the constraints are enforced. So  
// given values of a, b, c, d, and e at each of the points {i}, 
// this function constructs the affine transformation 
//
// phic = A*phi + B
//
// where A is a transformation matrix and B is the affine term. 
// If there are n points in the support of the shape function and m constraints, 
// there are (n-m) unconstrained points, and:
// 
// - A is an m x n matrix
// - phi is an n-dimensional column vector
// - phic and B are m-dimensional column vectors.
//
// Arguments:
// N - The polynomial least squares shape function.
// ghost_indices         - An array containing the indices of the ghost points
//                         on the shape function's present domain.
// num_ghosts            - The number of ghost points, and the length of ghost_indices.
//                         This number must be less than the number of points in the domain
//                         of the shape function.
// constraint_points     - An array containing the points in space at which the 
//                         constraints are enforced.
// a, b, c, d, e         - Arrays (of length num_constraints) whose ith entries contain
//                         the corresponding coefficients in the aforementioned constraint 
//                         equations.
// A                     - On output, this stores the matrix A in the affine 
//                         transformation above, in column-major order.
// B                     - On output, this stores the vector B in the affine 
//                         transformation above.
void poly_ls_shape_compute_ghost_transform(poly_ls_shape_t* N, int* constraint_indices, int num_constraints,
                                           point_t* constraint_points,
                                           double* a, double* b, double* c, double* d, double* e,
                                           double* A, double* B);

#endif
