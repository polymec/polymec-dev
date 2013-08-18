// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_GIFTWRAP_HULL_H
#define POLYMEC_GIFTWRAP_HULL_H

#include "core/polymec.h"


// Traverses the given (2D planar) points of a polygonal facet along their 
// convex hull, using the Giftwrap algorithm. Indices defining the ordering 
// are written to the indices array. The number of points that belong to the 
// convex hull is stored in count.
void giftwrap_hull(double* points, int num_points, int* indices, int* count);

// This version of giftwrap_hull also computes the area of the convex hull 
// using the fan algorithm.
void giftwrap_hull_with_area(double* points, int num_points, int* indices, int* count, double* area);

#endif

