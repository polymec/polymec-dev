// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "core/mesh.h"
#include "core/point_cloud.h"

// This type represents a mechanism for generating least-squares polynomial 
// fits for multi-component quantities on discrete domains.
typedef struct polynomial_fit_t polynomial_fit_t;

// This virtual table defines the behavior of a polynomial fit.
typedef struct 
{
  // Gets the number of interior neighbors associated with the point at the 
  // given index.
  int (*num_interior_neighbors)(void* context, int point_index);

  // Retrieves the indices of the interior neighbors of the point with the 
  // given index.
  void (*get_interior_neighbors)(void* context, int point_index, int* neighbor_indices);

  // Retrieves the coordinates and values of the given component at 
  // interior points in the vicinity of the point with the given index. The 
  // The coordinates of the points are stored in the array point_coords,
  // and the values of the component at each of the points are stored in the 
  // point_values array.
  void (*get_interior_data)(void* context, real_t* data, int component, 
                            int* point_indices, int num_points, 
                            point_t* point_coords, real_t* point_values);

  // Gets the number of boundary neighbors associated with the point at the 
  // given index.
  int (*num_boundary_neighbors)(void* context, int point_index);

  // Retrieves the indices of the boundary neighbors of the point with the 
  // given index.
  void (*get_boundary_neighbors)(void* context, int point_index, int* neighbor_indices);

  // Retrieves the coordinates and values of the given component at 
  // boundary points in the vicinity of the point with the given index. The 
  // The coordinates of the points are stored in the array point_coords,
  // the normal vectors on the boundary at each of the points are stored in 
  // the array boundary_normals. Note that point values are not retrieved,
  // as these are often determined by boundary conditions during the actual
  // fitting process (as represented by fit_component() below).
  void (*get_boundary_data)(void* context, real_t* data, int component, 
                            int* point_indices, int num_points, 
                            point_t* point_coords, vector_t* boundary_normals);

  // Returns the targeted degree of accuracy for the polynomial fit, given 
  // a number of points in the fit "stencil". This method effectively 
  // determines how aggressively a polynomial fit will pursue higher-order 
  // approximations to the solution.
  int (*targeted_degree)(void* context, int num_points);

  // Fits the data of a given component to a polynomial of the given degree, 
  // given interior and boundary data for points in the vicinity, storing the 
  // coefficients of the polynomial in the poly_coeffs array (in the order 
  // used by the polynomial_t type).
  void (*fit_component)(void* context, int component, int degree,
                        point_t* interior_points, real_t* interior_values, int num_interior_points,
                        point_t* boundary_points, vector_t* boundary_normals, int num_boundary_points,
                        real_t* poly_coeffs);

  // Destroys the context pointer.
  void (*dtor)(void* context);
} polynomial_fit_vtable;

// Constructs a generic polynomial fit that uses the provided functions to 
// compute fits on discrete domains.
polynomial_fit_t* polynomial_fit_new(const char* name,
                                     void* context,
                                     polynomial_fit_vtable vtable,
                                     int num_comps);

// Creates a polynomial fit that fits cell-centered data on a mesh at the 
// given degree. A fit_component() function implementation is required, as 
// well as a context it will be invoked with (and an associated destructor).
polynomial_fit_t* cc_fixed_degree_polynomial_fit_new(int num_comps,
                                                     mesh_t* mesh,
                                                     int degree,
                                                     void (*fit_component)(void*, int, int, point_t*, real_t*, int, point_t*, vector_t*, int, real_t*),
                                                     void* fit_component_context,
                                                     void (*dtor)(void*));

// Creates a polynomial fit that fits cell-centered data on a mesh at a 
// variable degree based on the number of neighbors it finds within a given 
// "depth." A fit_component() function is required.
polynomial_fit_t* cc_variable_degree_polynomial_fit_new(int num_comps,
                                                        mesh_t* mesh,
                                                        int neighbor_search_depth,
                                                        void (*fit_component)(void*, int, int, point_t*, real_t*, int, point_t*, vector_t*, int, real_t*),
                                                        void* fit_component_context,
                                                        void (*dtor)(void*));

// Creates a polynomial fit that fits point cloud data at the given degree.
// A fit_component() function is required.
polynomial_fit_t* point_cloud_fixed_degree_polynomial_fit_new(int num_comps,
                                                              point_cloud_t* points,
                                                              point_cloud_neighbor_search_t* search,
                                                              int degree,
                                                              void (*fit_component)(void*, int, int, point_t*, real_t*, int, point_t*, vector_t*, int, real_t*),
                                                              void* fit_component_context,
                                                              void (*dtor)(void*));

// Creates a polynomial fit that fits point cloud data at a variable degree 
// based on the number of neighbors it finds within a given "depth." A 
// fit_component() function is required.
polynomial_fit_t* point_cloud_variable_degree_polynomial_fit_new(int num_comps,
                                                                 point_cloud_t* points,
                                                                 point_cloud_neighbor_search_t* search,
                                                                 int neighbor_search_depth,
                                                                 void (*fit_component)(void*, int, int, point_t*, real_t*, int, point_t*, vector_t*, int, real_t*),
                                                                 void* fit_component_context,
                                                                 void (*dtor)(void*));

// Destroys the given polynomial fit.
void polynomial_fit_free(polynomial_fit_t* fit);

// Returns the number of components in the quantity represented by the 
// polynomial fit.
int polynomial_fit_num_comps(polynomial_fit_t* fit);

// Computes the polynomial fit at the given indexed point in the discrete 
// domain. This discrete domain can be a mesh, a point_cloud, or another 
// abstract discrete domain, as long as neighbor relationships are defined 
// for points.
void polynomial_fit_compute(polynomial_fit_t* fit, real_t* data, int point_index);

// Returns the degree of the present polynomial fit. Returns -1 if no fit 
// has yet been computed.
int polynomial_fit_degree(polynomial_fit_t* fit);

// Evaluate the computed polynomial fit at the given point in space.
void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value);

// Evaluate the given partial derivative of the computed polynomial fit at 
// the given point in space.
void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               point_t* x, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               real_t* deriv);

#endif
