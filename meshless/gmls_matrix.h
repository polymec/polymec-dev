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
#include "model/stencil.h"
#include "model/point_weight_function.h"
#include "meshless/gmls_functional.h"

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

  // This function computes the vector y corresponding to the effective 
  // displacement from a point xi to a point xi. This vector is to be handed 
  // to a point_weight_function to determine the MLS weight for point xj as 
  // it contributes to a MLS approximation on the subdomain i (centered about
  // xi). If this method is not provided, y = xj - xi.
  void (*compute_weight_displacement)(void* context, 
                                      int i, point_t* xi, 
                                      int j, point_t* xj, 
                                      vector_t* y);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); // Destructor
} gmls_matrix_vtable;

// Creates a generalized MLS functional with the given name, context, virtual 
// table, (multicomponent) polynomial basis, and weight function. The 
// weight function is consumed by the matrix. The solution is assumed 
// to be a vector of degrees of freedom with the given number of components, 
// expressed in node-major order.
gmls_matrix_t* gmls_matrix_new(const char* name,
                               void* context,
                               gmls_matrix_vtable vtable,
                               multicomp_poly_basis_t* poly_basis,
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

// Computes the given row of the GMLS matrix in order to enforce a (collocated) 
// Dirichlet boundary condition for the solution there. Column indices and 
// coefficients are stored in the columns and coeffs arrays.
void gmls_matrix_compute_dirichlet_row(gmls_matrix_t* matrix,
                                       int row,
                                       int* columns,
                                       real_t* coeffs);

// Computes the given row of the GMLS matrix in order to enforce a (collocated) 
// Neumann boundary condition for the solution there. The given operator op, 
// represented as a space-time function, is applied as a transformation to the 
// gradient of the polynomial basis at the boundary. It can have 3 or 6 
// components: 
// * If it has 3, it is interpreted as a vector to be dotted into the 
//   polynomial gradient.
// * If it has 6, it is interpreted as a symmetric tensor to be dotted into the 
//   polynomial gradient.
// * If it has 9, it is interpreted as a nonsymmetric tensor to be dotted into 
//   the polynomial gradient.
// The form of the boundary condition reads op(grad u) = value.
// Column indices and coefficients are stored in the columns and coeffs arrays.
void gmls_matrix_compute_neumann_row(gmls_matrix_t* matrix,
                                     int row,
                                     st_func_t* op,
                                     int* columns,
                                     real_t* coeffs);

// Computes the given row of the GMLS matrix in order to enforce a (collocated) 
// Robin boundary condition for the solution there. The given operators op1 and 
// op2, represented as a space-time functions, are applied as transformations 
// to the value and gradient of the polynomial basis at the boundary. op1, 
// which multiplies the value of the polynomial basis itself, must have exactly 
// one component. op2 can have 3 or 6 components: 
// * If it has 3, it is interpreted as a vector to be dotted into the 
//   polynomial gradient.
// * If it has 6, it is interpreted as a symmetric tensor to be dotted into the 
//   polynomial gradient.
// The form of the boundary condition reads op1(u) + op2(grad u) = value.
// Column indices and coefficients are stored in the columns and coeffs arrays.
void gmls_matrix_compute_robin_row(gmls_matrix_t* matrix,
                                   int row,
                                   st_func_t* alpha_op,
                                   st_func_t* beta_op,
                                   int* columns,
                                   real_t* coeffs);

// Creates a GMLS matrix using a point cloud and a given stencil to provide information 
// about nodes contributing to subdomains. The point cloud and the stencil are both 
// borrowed by the matrix, not consumed. The weight function displacement 
// y = |x1 - x2|/R1, where R1 is the extent of the subdomain corresponding 
// to x1.
gmls_matrix_t* stencil_based_gmls_matrix_new(multicomp_poly_basis_t* poly_basis,
                                             point_weight_function_t* W,
                                             point_cloud_t* points,
                                             real_t* extents,
                                             stencil_t* stencil);

#endif
