// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POINT_CLOUD_FIELD_H
#define POLYMEC_POINT_CLOUD_FIELD_H

#include "core/point_cloud.h"

typedef struct point_cloud_field_t point_cloud_field_t;

// Constructs a new point cloud field with the given number of components
// on the given point cloud.
point_cloud_field_t* point_cloud_field_new(point_cloud_t* cloud,
                                           int num_components);

// Destroys the given point cloud field.
void point_cloud_field_free(point_cloud_field_t* field);

// Returns the point cloud on which the field is defined.
point_cloud_t* point_cloud_field_cloud(point_cloud_field_t* field);

// Returns the number of components in the field.
int point_cloud_field_num_components(point_cloud_field_t* field);

// Returns the number of local values in the field.
int point_cloud_field_num_local_values(point_cloud_field_t* field);

// Returns the number of ghost values in the field.
int point_cloud_field_num_ghost_values(point_cloud_field_t* field);

// Returns the field's data.
real_t* point_floud_field_data(point_cloud_field_t* field);

// Defines a 2D array that allows the field to be indexed thus:
// f[i][c] returns the cth component of the value for the ith point.
#define DECLARE_POINT_CLOUD_FIELD_ARRAY(array, field) \
DECLARE_2D_ARRAY(real_t, array, point_cloud_field_data(field), point_cloud_field_num_local_values(field) + point_cloud_field_num_ghost_values(field), point_cloud_field_num_components(field))

#endif

