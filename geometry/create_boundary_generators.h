// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_BOUNDARY_GENERATORS_H
#define POLYMEC_CREATE_BOUNDARY_GENERATORS_H

#include "core/point.h"
#include "core/array.h"

// This function creates a set of generators that can be used as stationary 
// generators in a Centroidal Voronoi Tessellation (CVT) algorithm.
void create_boundary_generators(ptr_array_t* surface_points, 
                                ptr_array_t* surface_normals, 
                                ptr_array_t* surface_tags,
                                point_t** boundary_generators,
                                int* num_boundary_generators,
                                char*** tag_names,
                                int_array_t*** tags,
                                int* num_tags);

#endif

