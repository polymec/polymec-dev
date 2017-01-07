// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_BASIS_H
#define POLYMEC_POINT_BASIS_H

#include "core/point.h"

// This class represents a set of basis functions defined on a point cloud
// within a neighborhood. These basis functions can represent quantities on 
// the point cloud.
typedef struct point_basis_t point_basis_t;

// Here's a virtual table used to define the behavior of a point function.
typedef struct
{
  // Returns the number of points in the neighborhood of the point i.
  int (*neighborhood_size)(void* context, int i);
  // This method gets the points in the neighborhood of the point i.
  void (*get_neighborhood_points)(void* context, int i, point_t* points);
  // This (optional) method does any work associated with precomputing data
  // for evaluating the shape function in the neighborhood of point i.
  void (*set_neighborhood)(void* context, int i);
  // This method computes the values of the shape functions on each of the 
  // points in the neighborhood of the point i, evaluated at the point x. If 
  // the gradient argument is non-NULL, the gradient of the shape function is 
  // also computed at these points.
  void (*compute)(void* context, int i, point_t* x, real_t* values, vector_t* gradients);
  // This destructor destroys the context.
  void (*dtor)(void* context);
} point_basis_vtable;

// Creates an instance of a shape function with the given name and 
// characteristics. 
point_basis_t* point_basis_new(const char* name, 
                               void* context, 
                               point_basis_vtable vtable);

// Destroys the shape function.
void point_basis_free(point_basis_t* phi);

// Sets the point within the domain in whose vicinity the shape function 
// will be defined, using the stencil for that point to define the neighborhood.
void point_basis_set_neighborhood(point_basis_t* phi, int point_index);

// Returns the number of points in the current neighborhood, or -1 if the 
// neighborhood has not been set with point_basis_set_neighborhood.
int point_basis_num_points(point_basis_t* phi);

// Fetches the points in the current neighborhood, filling the points array
// (of length point_basis_num_points(phi)).
void point_basis_get_points(point_basis_t* phi, point_t* points);

// Computes the values of the shape functions for the points in the current 
// neighborhood, evaluating them at the point x, filling the array with 
// values. If gradients is non-NULL, the gradients will also be computed.
void point_basis_compute(point_basis_t* phi, 
                         point_t* x,
                         real_t* values,
                         vector_t* gradients);

#endif
