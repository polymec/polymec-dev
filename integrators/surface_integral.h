// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SURFACE_INTEGRAL_H
#define SURFACE_INTEGRAL_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/st_func.h"

// This class represents an interface for approximating surface integrals 
// using a given quadrature rule / surface. Objects of this type are 
// garbage-collected.
typedef struct surface_integral_t surface_integral_t;

// The following virtual table allows the implementation of a surface integral.
typedef struct
{
  // This function returns the number of quadrature points for the quadrature 
  // rule on the ith subdomain in a spatial discretization.
  int (*num_quad_points)(void* context, int i);

  // This function computes the quadrature points, weights, and normals for 
  // the ith subdomain in a spatial discretization.
  void (*get_quadrature)(void* context, int i, point_t* points, real_t* weights, vector_t* normals);

  // Destructor function.
  void (*dtor)(void* context);
} surface_integral_vtable;

// Construct a new surface integral with the given name, context, and 
// virtual table.
surface_integral_t* surface_integral_new(const char* name,
                                         void* context,
                                         surface_integral_vtable vtable);

// Sets up the surface integral to compute on the ith subdomain in a spatial 
// discretization. This could be a cell within a mesh, or a point within a 
// point cloud. This must be called before any integrals are computed.
void surface_integral_set_domain(surface_integral_t* integ, int i);

// Computes the surface integral of the given spatial function f, storing it 
// in integral, which should be equal in size to the number of components in f.
void surface_integral_compute(surface_integral_t* integ, sp_func_t* f, real_t* integral);

// Computes the surface integral of the given space-time function f at time t, 
// storing it in integral, which should be equal in size to the number of 
// components in f.
void surface_integral_compute_at_time(surface_integral_t* integ, real_t t, st_func_t* f, real_t* integral);

// Returns the surface integral of the given 3D vector-valued spatial function 
// f dotted into the normal vector of the surface.
real_t surface_integral_dot(surface_integral_t* integ, sp_func_t* f);

// Returns the surface integral of the given 3D vector-valued space-time 
// function f dotted into the normal vector of the surface at time t.
real_t surface_integral_dot_at_time(surface_integral_t* integ, 
                                    real_t t, st_func_t* f);

// Returns the number of quadrature points in the underlying approximation 
// of the surface integral.
int surface_integral_num_points(surface_integral_t* integ);

// Retrieves the quadrature points, weights, and normal vectors in the 
// underlying approximation, storing them in the points, weights, and normals 
// arrays, which should be large enough to store them.
void surface_integral_get_quadrature(surface_integral_t* integ,
                                     point_t* points,
                                     real_t* weights,
                                     vector_t* normals);

// Traverses the quadrature points in the surface integral, storing the 
// point, weight, and normal vector in the given pointers. Returns true if 
// the traversal yielded another quadrature point/weight/normal, false if it 
// is finished. Set *pos = 0 to reset the traversal.
bool surface_integral_next_quad_point(surface_integral_t* integ, int* pos,
                                      point_t* point, real_t* weight, vector_t* normal);

#endif

