// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_CONVEX_HULL_H
#define POLYMEC_CREATE_CONVEX_HULL_H

#include "core/point.h"

// Creates the convex hull of the given points, returning it in an array and 
// setting *hull_size to the number of points in it.
point_t* create_convex_hull(point_t* points, int num_points, int* hull_size);
 
#endif

