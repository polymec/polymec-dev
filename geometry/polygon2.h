// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYGON2_H
#define POLYMEC_POLYGON2_H

#include "core/point2.h"

// This class represents a polygon in the plane. Objects of this type are 
// garbage-collected.
typedef struct polygon2_t polygon2_t;

// Creates a new polygon in the plane given an ordered sequence of points 
// that identify consecutive vertices.
polygon2_t* polygon2_new(point2_t* vertices, int num_vertices);

// Creates a new polygon in the plane given a set of vertices and an 
// ordering in which they are traversed.
polygon2_t* polygon2_new_with_ordering(point2_t* points, int* ordering, int num_points);

// Creates a new convex polygon in the plane by applying the "gift-wrap" 
// convex hull algorithm to the set of points.
polygon2_t* polygon2_giftwrap(point2_t* points, int num_points);

// Creates a new nonconvex "star-shaped" polygon in the plane by ordering 
// the given set of points according to their angles of displacement about 
// a given center x0. The points must all be distinct from x0.
polygon2_t* polygon2_star(point2_t* x0, point2_t* points, int num_points);

// Returns the number of vertices in the polygon.
int polygon2_num_vertices(polygon2_t* poly);

// Returns the ordering of the vertices in the polygon in terms of 
// the original array of vertices given in the constructor.
int* polygon2_ordering(polygon2_t* poly);

// Allows the traversal of the vertices in the polygon.
bool polygon2_next_vertex(polygon2_t* poly, int* pos, point2_t** vertex);

// Returns the area of the polygon.
real_t polygon2_area(polygon2_t* poly);

// Clones the polygon, returning an exact copy.
polygon2_t* polygon2_clone(polygon2_t* poly);

// Clips the given polygon by intersecting it with another. Here, poly is 
// the polygon that is to be clipped, and it is modified in place.
// This algorithm only works reliably if both polygons are convex.
void polygon2_clip(polygon2_t* poly, polygon2_t* other);

#endif

