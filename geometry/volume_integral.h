// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef VOLUME_INTEGRAL_H
#define VOLUME_INTEGRAL_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/st_func.h"

// This class represents an interface for approximating volume integrals 
// using a given quadrature rule / volume. Objects of this type are 
// garbage-collected.
typedef struct volume_integral_t volume_integral_t;

// The following virtual table allows the implementation of a volume integral.
typedef struct
{
  // This function returns the number of quadrature points for the quadrature 
  // rule on the ith subdomain in a spatial discretization.
  int (*num_quad_points)(void* context, int i);

  // This function computes the quadrature points and weights for the ith 
  // subdomain in a spatial discretization.
  void (*get_quadrature)(void* context, int i, point_t* points, real_t* weights);

  // Destructor function.
  void (*dtor)(void* context);
} volume_integral_vtable;

// Construct a new volume integral with the given name, context, and 
// virtual table.
volume_integral_t* volume_integral_new(const char* name,
                                       void* context,
                                       volume_integral_vtable vtable);

// Sets up the volume integral to compute on the ith subdomain in a spatial 
// discretization. This could be a cell within a mesh, or a point within a 
// point cloud. This must be called before any integrals are computed.
void volume_integral_set_domain(volume_integral_t* integ, int i);

// Returns the volume integral of the given spatial function f, storing it 
// in integral, which should be equal in size to the number of components in f.
void volume_integral_compute(volume_integral_t* integ, sp_func_t* f, real_t* integral);

// Returns the volume integral of the given space-time function f at time t, 
// storing it in integral, which should be equal in size to the number of 
// components in f.
void volume_integral_compute_at_time(volume_integral_t* integ, real_t t, st_func_t* f, real_t* integral);

// Returns the number of quadrature points in the underlying approximation 
// of the volume integral.
int volume_integral_num_points(volume_integral_t* integ);

// Retrieves the quadrature points and weights in the underlying approximation,
// storing them in the points and weights arrays, which should be large 
// enough to store them.
void volume_integral_get_quadrature(volume_integral_t* integ,
                                    point_t* points,
                                    real_t* weights);

// Traverses the quadrature points in the volume integral, storing the 
// point and weight in the given pointers. Returns true if the traversal 
// yielded another quadrature point/weight, false if it is finished. Set 
// *pos = 0 to reset the traversal.
bool volume_integral_next_quad_point(volume_integral_t* integ, int* pos,
                                     point_t* point, real_t* weight);

#endif

