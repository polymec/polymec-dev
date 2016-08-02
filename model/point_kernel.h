// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_KERNEL_H
#define POLYMEC_POINT_KERNEL_H

#include "core/point.h"

// Here's a "kernel" that can be used to construct approximations of quantities
// on point clouds. By "kernel" we mean a distribution function W(x, y, h) that 
// becomes a "delta function" as the value of the "kernel size" h approaches 
// zero. Objects of this type are garbage-collected.
typedef struct point_kernel_t point_kernel_t;

// Creates a new point kernel with the behavior defined by the context, 
// the compute function, and a destructor.
point_kernel_t* point_kernel_new(const char* name,
                                 void* context,
                                 void (*compute)(void* context, point_t* points, real_t* kernel_sizes, int num_points, point_t* x, real_t* values, vector_t* gradients),
                                 void (*dtor)(void* context));

// Evaluates the kernels centered on the given points (with the given 
// "kernel sizes"), computing their values and (if gradients != NULL) their 
// gradients at the point x.
void point_kernel_compute(point_kernel_t* kernel, 
                          point_t* points,
                          real_t* kernel_sizes,
                          int num_points, 
                          point_t* x, 
                          real_t* values,
                          vector_t* gradients);

// This is a simple cubic B-spline kernel of the form
// W(x, y, h) = 1 + 6 * q**2 * (q-1), 0 <= q < 1/2,
//            = 2 * (1-q)**3,         1/2 <= q < 1,
//            = 0,                    q >= 1 
// where q = ||x-y||/h. This kernel is commonly used in SPH.
point_kernel_t* cubic_bspline_point_kernel_new();
              
// This is a simple fourth-order spline-based kernel of the form
// W(x, y, h) = 1 - 6*(q/q_max)**2 + 8*(q/q_max)**3 - 3*(q/q_max)**4 
// for 0 <= q <= q_max, 0 for q > q_max. Here, q = ||x-y||/h.
point_kernel_t* spline4_point_kernel_new(real_t q_max);
              
#endif
