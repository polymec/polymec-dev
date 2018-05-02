// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_CLOUD_FIELD_H
#define POLYMEC_POINT_CLOUD_FIELD_H

#include "geometry/point_cloud.h"

// This thin wrapper represents a field of values defined on the points of 
// a point cloud.
typedef struct 
{
  // The underlying point cloud.
  point_cloud_t* cloud;

  // The number of components for a datum in the field.
  size_t num_components;

  // The number of locally-stored values in the field.
  size_t num_local_values;

  // The number of ghost values in the field (which changes with the 
  // number of ghost points in the cloud).
  size_t num_ghost_values;

  // Data for the field, and its storage capacity.
  real_t* data;
  size_t capacity;
} point_cloud_field_t;

// Constructs a new point cloud field with the given number of components
// on the given point cloud.
point_cloud_field_t* point_cloud_field_new(point_cloud_t* cloud,
                                           size_t num_components);

// Destroys the given point cloud field.
void point_cloud_field_free(point_cloud_field_t* field);

// Defines a 2D array that allows the field to be indexed thus:
// f[i][c] returns the cth component of the value for the ith point.
#define DECLARE_POINT_CLOUD_FIELD_ARRAY(array, field) \
DECLARE_2D_ARRAY(real_t, array, field->data, field->num_local_values + field->num_ghost_values, field->num_components)

#endif

