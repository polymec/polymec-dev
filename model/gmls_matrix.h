// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GMLS_MATRIX_H
#define POLYMEC_GMLS_MATRIX_H

#include "core/polynomial.h"
#include "core/point_cloud.h"
#include "model/gmls_functional.h"
#include "model/point_weight_function.h"
#include "model/stencil.h"

// This class uses a Generalized Moving Least Squares (GMLS) functional to construct 
// a matrix such that 
//
// u(x, t) = A(x, t) * U(t).
//
// where A is a matrix that maps the nodal values of u (stored in the vector U) to its 
// value at x and t.
// 
// This construct is useful in constructing Meshless Local Petrov-Galerkin 
// (MLPG) methods that do not rely on shape functions, such as the Direct 
// MLPG (DMLPG) method.
typedef struct gmls_matrix_t gmls_matrix_t;

// This virtual table defines the behavior of a GMLS matrix, essentially telling it 
// how to get nodes that contribute to a given subdomain. We provide this flexibility 
// to allow for dynamically changing connectivity.
typedef struct
{
  // This function returns the number of nodes contributing to the ith 
  // subdomain.
  int (*num_nodes)(void* context, int i);

  // This function retrieves the indices of the nodes contributing to the 
  // ith subdomain, placing them into the nodes array.
  void (*get_nodes)(void* context, int i, int* nodes);

  // This function retrieves the locations of the nodes with the given indices.
  void (*get_points)(void* context, int* nodes, int num_nodes, point_t* points);

  // This function retrieves the distance scale factors associated with the 
  // nodes with the given indices. These scale factors are used to scale the 
  // weights associated with the nodes. If this method is not provided, 
  // the weight function will scale the distances between points.
  void (*get_scale_factors)(void* context, int* nodes, int num_nodes, real_t* scale_factors);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); // Destructor
} gmls_matrix_vtable;

// Creates a generalized MLS functional with the given name, context, virtual 
// table, (multicomponent) polynomial basis, and weight function. The 
// weight function is consumed by the matrix. The solution is assumed 
// to be a vector of degrees of freedom expressed in node-major order.
gmls_matrix_t* gmls_matrix_new(const char* name,
                               void* context,
                               gmls_matrix_vtable vtable,
                               point_weight_function_t* W);

// Destroys the given GMLS matrix.
void gmls_matrix_free(gmls_matrix_t* matrix);

// Returns the number of nonzero columns for the given row in the matrix.
int gmls_matrix_num_columns(gmls_matrix_t* matrix, int row);

// Evaluates the given row of the GMLS matrix at using the given functional 
// evaluated at time t, placing the indices of the nonzero columns into the 
// columns array, and the coefficients into the coefficients array. These 
// arrays should both sized using gmls_matrix_num_columns(matrix, row).
void gmls_matrix_compute_row(gmls_matrix_t* matrix,
                             int row,
                             gmls_functional_t* lambda,
                             real_t t,
                             int* columns,
                             real_t* coeffs);

// Creates a GMLS matrix using a point cloud and a given stencil to provide information 
// about nodes contributing to subdomains. The point cloud and the stencil are both 
// borrowed by the matrix, not consumed.
gmls_matrix_t* stencil_based_gmls_matrix_new(point_weight_function_t* W,
                                             point_cloud_t* points,
                                             stencil_t* stencil);

#endif