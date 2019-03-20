// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LEAST_SQUARES_H
#define POLYMEC_LEAST_SQUARES_H

#include "core/polymec.h"
#include "core/point.h"

// These functions and types are intended to help with the construction of
// least-squares fits to scattered data.

/// \addtogroup core core
///@{

/// Computes the coefficients A and b for a linear regression y = A*x + B
/// given arrays of (scalar) x and y values. Also computes the variance sigma.
void linear_regression(real_t* x, real_t* y, int N, real_t* A, real_t* B, real_t* sigma);

/// \class ls_weight_func
/// This type is a weight function used for least-squares systems.
/// \refcounted
typedef struct ls_weight_func_t ls_weight_func_t;

/// This is the signature for least-squares weight functions. Arguments:
/// - A context pointer for this weighting function.
/// - The displacement y = x - x0 of the evaluated point from the function's "center".
/// - A pointer to storage for the weighting function value.
/// - A pointer to storage for the weighting function gradient.
typedef void (*ls_weight_func_eval_func)(void* context, vector_t* y, real_t* W, vector_t* grad_W);

/// This is the signature for setting the "domain" (set of points) that will
/// define the extent of a least-squares weight function.
typedef void (*ls_weight_func_set_domain_func)(void* context, point_t* x0, point_t* points, int num_points);

/// Destructor function for weight function context.
typedef void (*ls_weight_func_dtor)(void*);

/// \struct ls_weight_func_vtable
/// This virtual table must be implemented for a least-squares weight function.
typedef struct
{
  ls_weight_func_set_domain_func set_domain;
  ls_weight_func_eval_func       eval;
  ls_weight_func_dtor            dtor;
} ls_weight_func_vtable;

/// Creates a least-squares weight function from a name, a context, and a
/// virtual table.
/// \memberof ls_weight_func
ls_weight_func_t* ls_weight_func_new(const char* name,
                                     void* context,
                                     ls_weight_func_vtable vtable);

/// Returns the name of the weight function.
/// \memberof ls_weight_func
const char* ls_weight_func_name(ls_weight_func_t* W);

/// Sets the domain (the center point x0 and the points on which data exist) of the weight function.
/// \memberof ls_weight_func
void ls_weight_func_set_domain(ls_weight_func_t* W, point_t* x0, point_t* points, int num_points);

/// Evaluates the value and the gradient of the weight function at the given point x.
/// \memberof ls_weight_func
void ls_weight_func_eval(ls_weight_func_t* W, point_t* x, real_t* value, vector_t* gradient);

/// Returns the current center point x0 of the weight function.
/// \memberof ls_weight_func
point_t* ls_weight_func_x0(ls_weight_func_t* W);

/// Computes the least squares system for a pth-order polynomial fit to a set
/// of scattered point data, centered about the point x0. This computes the
/// moment matrix and the right-hand side vector for a polynomial fit to the
/// given (scalar) data on the given points.
/// \memberof ls_weight_func
void compute_poly_ls_system(int p, point_t* x0, point_t* points, int num_points,
                            real_t* data, real_t* moment_matrix, real_t* rhs);

/// This computes the weighted least squares system for the weighting function W, centered
/// about the point x0 and with spatial extent h.
/// \memberof ls_weight_func
void compute_weighted_poly_ls_system(int p, ls_weight_func_t* W, point_t* x0,
                                     point_t* points, int num_points, real_t* data,
                                     real_t* moment_matrix, real_t* rhs);

///@}

#endif
