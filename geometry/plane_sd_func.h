// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PLANE_SD_FUNC_H
#define POLYMEC_PLANE_SD_FUNC_H

#include "geometry/sd_func.h"

// This signed distance function represents a plane with a given normal 
// vector and point.
sd_func_t* plane_sd_func_new(vector_t* n, point_t* x);

// Construct a plane that contains the three points p1, p2, and p3, 
// assuming that these points are not co-linear.
sd_func_t* plane_sd_func_from_points(point_t* p1, point_t* p2, point_t* p3);

#endif

