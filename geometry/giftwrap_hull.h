#ifndef POLYMEC_GIFTWRAP_HULL_H
#define POLYMEC_GIFTWRAP_HULL_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif


// Traverses the given (2D planar) points of a polygonal facet along their 
// convex hull, using the Giftwrap algorithm. Indices defining the ordering 
// are written to the indices array. The number of points that belong to the 
// convex hull is stored in count.
void giftwrap_hull(double* points, int num_points, int* indices, int* count);

// This version of giftwrap_hull also computes the area of the convex hull 
// using the fan algorithm.
void giftwrap_hull_with_area(double* points, int num_points, int* indices, int* count, double* area);

#ifdef __cplusplus
}
#endif

#endif

