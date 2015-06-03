// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SHAPE_FUNCTION_H
#define POLYMEC_SHAPE_FUNCTION_H

#include "core/point.h"

// Here's a "kernel" that can be used to construct shape functions.
// This interface is not strictly necessary for defining a given shape 
// function, but is provided as a convenience here.
typedef struct shape_function_kernel_t shape_function_kernel_t;

// Creates a kernel function that can be used to construct a shape function.
shape_function_kernel_t* shape_function_kernel_new(const char* name,
                                                   void* context,
                                                   void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* values, vector_t* gradients),
                                                   void (*dtor)(void* context));

// Destroys the given shape function kernel.
void shape_function_kernel_free(shape_function_kernel_t* kernel);

// Evaluates the kernel functions centered on the given points (with "extents"), computing 
// their values and (if gradients != NULL) their gradients at the point x.
void shape_function_kernel_compute(shape_function_kernel_t* kernel, 
                                   point_t* points,
                                   real_t* extents,
                                   int num_points, 
                                   point_t* x, 
                                   real_t* values,
                                   vector_t* gradients);
              
// This class represents a shape function defined on a set of points in a 
// neighborhood within a point cloud.
typedef struct shape_function_t shape_function_t;

// Here's a virtual table used to define the behavior of a shape function.
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
} shape_function_vtable;

// Creates an instance of a shape function with the given name and 
// characteristics. 
shape_function_t* shape_function_new(const char* name, 
                                     void* context, 
                                     shape_function_vtable vtable);

// Destroys the shape function.
void shape_function_free(shape_function_t* phi);

// Sets the point within the domain in whose vicinity the shape function 
// will be defined, using the stencil for that point to define the neighborhood.
void shape_function_set_neighborhood(shape_function_t* phi, int point_index);

// Returns the number of points in the current neighborhood, or -1 if the 
// neighborhood has not been set with shape_function_set_neighborhood.
int shape_function_num_points(shape_function_t* phi);

// Fetchs the points in the current neighborhood, filling the points array
// (of length shape_function_num_points(phi)).
void shape_function_get_points(shape_function_t* phi, point_t* points);

// Computes the values of the shape functions for the points in the current 
// neighborhood, evaluating them at the point x, filling the array with 
// values. If gradients is non-NULL, the gradients will also be computed.
void shape_function_compute(shape_function_t* phi, 
                            point_t* x,
                            real_t* values,
                            vector_t* gradients);

#endif
