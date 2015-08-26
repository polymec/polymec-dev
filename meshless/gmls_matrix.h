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

// This class represents an operator that is applied to a solution at the boundary 
// to represent part or all of an expression for a boundary condition. It may be a 
// scalar or tensor operator, time dependent or independent, and nonlinear (dependent 
// on the solution itself) or linear (not so). Objects of this type are garbage-collected.
typedef struct gmls_matrix_bc_operator_t gmls_matrix_bc_operator_t;

// This virtual table defines an interface that describes the behavior of a boundary 
// condition operator.
typedef struct
{
  // This function computes the value of the operator at (x, t), and (for nonlinear 
  // operators only) given the value of the solution. It computes a 
  void (*compute)(void* context, point_t* x, real_t t, real_t* solution, real_t* operator_matrix);

  // This is a destructor that destroys the given context.
  void (*dtor)(void* context); 
} gmls_matrix_bc_operator_vtable;

// Creates a boundary condition operator, given a name, a context, a virtual table,
// and corresponding to the given number of solution components.
gmls_matrix_bc_operator_t* gmls_bc_operator_new(const char* name, 
                                                void* context,
                                                gmls_matrix_bc_operator_vtable vtable,
                                                int num_components);

// Returns the number of components for the solutions compatible with this operator.
int gmls_matrix_bc_operator_num_components(gmls_matrix_bc_operator_t* op);

// Computes an NxN operator matrix to be used in the application of boundary conditions,
// evaluated at the point x and the time t and (for nonlinear operators) the given 
// solution. The resulting matrix is stored in column-major format and placed into the 
// op_matrix array.
void gmls_matrix_bc_operator_compute(gmls_matrix_bc_operator_t* op, 
                                     point_t* x,
                                     real_t t,
                                     real_t* solution,
                                     real_t* op_matrix);

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

// Returns the number of matrix coefficients that correspond to the ith 
// node in the GMLS approximation.
int gmls_matrix_num_coeffs(gmls_matrix_t* matrix, int i);

// Evaluates the coefficients of the GMLS matrix for the solution components 
// at node i, using the given functional evaluated at time t. The indices 
// of the rows and columns for these coefficients are placed into the given 
// rows and columns arrays, and the coefficients into the coefficients array. These 
// arrays should both sized using gmls_matrix_num_columns(matrix, row).
// The values of the solution may be given in a component-minor-ordered array
// for nonlinear problems; for linear problems solution may be set to NULL.
void gmls_matrix_compute_coeffs(gmls_matrix_t* matrix,
                                int i,
                                gmls_functional_t* lambda,
                                real_t t,
                                real_t* solution,
                                int* rows,
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

// Creates a quadrature rule that can be used to enforce boundary conditions
// for this matrix.
volume_integral_t* gmls_matrix_bc_quadrature_new(gmls_matrix_t* matrix);

// Creates a functional that can be used with this matrix to enforce a very 
// simple Dirichlet boundary condition on each component of the solution.
// NOTE: This functional only works properly on the standard multicomponent 
// NOTE: polynomial basis.
gmls_functional_t* gmls_matrix_dirichlet_bc_new(gmls_matrix_t* matrix);

// Creates a functional that can be used with this matrix to enforce a very 
// simple Neumann boundary condition on each component of the solution, 
// given a (3-component) boundary normal vector n(x, t).
// NOTE: This functional only works properly on the standard multicomponent 
// NOTE: polynomial basis.
gmls_functional_t* gmls_matrix_neumann_bc_new(gmls_matrix_t* matrix,
                                              st_func_t* n);

// Creates a functional that can be used with this matrix to enforce a very 
// simple Robin boundary condition on each component of the solution, 
// given a (3-component) boundary normal vector n(x, t) and coefficients 
// alpha and beta that multiply the solution components and their gradients, 
// respectively.
// NOTE: This functional only works properly on the standard multicomponent 
// NOTE: polynomial basis.
gmls_functional_t* gmls_matrix_robin_bc_new(gmls_matrix_t* matrix,
                                            st_func_t* n,
                                            real_t alpha,
                                            real_t beta);

#endif
