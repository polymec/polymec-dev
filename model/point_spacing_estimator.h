// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_SPACING_ESTIMATOR_H
#define POLYMEC_POINT_SPACING_ESTIMATOR_H

#include "core/point.h"
#include "model/stencil.h"

// This class provides a mechanism for estimating a "spacing" dx in the 
// vicinity of a given point within some discretization. Objects of this 
// type are garbage-collected.
typedef struct point_spacing_estimator_t point_spacing_estimator_t;

// This virtual table must be implemented by a point spacing estimator.
typedef struct 
{
  // Returns a grid spacing dx for the ith point.
  real_t (*dx)(void* context, int i);
  // Destructor.
  void (*dtor)(void* context);
} point_spacing_estimator_vtable;

// Constructs a new point cloud spacing estimator with the given name, 
// context, and virtual table, and associated with the given point cloud.
point_spacing_estimator_t* point_spacing_estimator_new(const char* name,
                                                       void* context,
                                                       point_spacing_estimator_vtable vtable);

// Returns the point spacing in the vicinity of the ith point for this 
// estimator.
real_t point_spacing_estimator_dx(point_spacing_estimator_t* estimator, int i);

// Shake-n-bake stencil-based point-spacing estimator. Given a point cloud
// a stencil object identifying neighbor relations between points, this estimates
// the average distance between the neighbors of a point. The point cloud and 
// stencil are not consumed by the estimator.
point_spacing_estimator_t* stencil_point_spacing_estimator_new(point_cloud_t* points,
                                                               stencil_t* stencil);

#endif

