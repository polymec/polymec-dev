#ifndef POLYMEC_POLYGON_H
#define POLYMEC_POLYGON_H

#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class represents a polygon embedded in 3-dimensional space. 
// Objects of this type are garbage-collected.
typedef struct polygon_t polygon_t;

// Creates a new polygon given an ordered sequence of points 
// that identify consecutive vertices.
polygon_t* polygon_new(point_t* vertices, int num_vertices);

// Creates a new polygon given a set of vertices and an ordering in which 
// they are traversed.
polygon_t* polygon_new_with_ordering(point_t* points, int* ordering, int num_points);

// Creates a new convex polygon by applying the "gift-wrap" 
// convex hull algorithm to the set of points.
polygon_t* polygon_giftwrap(point_t* points, int num_points);

// Returns the number of vertices in the polygon.
int polygon_num_vertices(polygon_t* poly);

// Returns the ordering of the vertices in the polygon in terms of 
// the original array of vertices given in the constructor.
int* polygon_ordering(polygon_t* poly);

// Allows the traversal of the vertices in the polygon.
bool polygon_next_vertex(polygon_t* poly, int* pos, point_t** vertex);

// Returns the area of the polygon.
double polygon_area(polygon_t* poly);

// Clones the polygon, returning an exact copy.
polygon_t* polygon_clone(polygon_t* poly);

// Clips the given polygon by intersecting it with another. Here, poly is 
// the polygon that is to be clipped, and it is modified in place.
// This algorithm only works reliably if both polygons are convex.
void polygon_clip(polygon_t* poly, polygon_t* other);

#ifdef __cplusplus
}
#endif

#endif

