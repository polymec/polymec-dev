// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GIFTWRAP_HULL_H
#define POLYMEC_GIFTWRAP_HULL_H

#include "core/polymec.h"


// Traverses the given (2D planar) points of a polygonal facet along their 
// convex hull, using the Giftwrap algorithm. Indices defining the ordering 
// are written to the indices array. The number of points that belong to the 
// convex hull is stored in count.
void giftwrap_hull(real_t* points, int num_points, int* indices, int* count);

// This version of giftwrap_hull also computes the area of the convex hull 
// using the fan algorithm.
void giftwrap_hull_with_area(real_t* points, int num_points, int* indices, int* count, real_t* area);

#endif

