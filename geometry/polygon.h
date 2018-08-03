// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYGON_H
#define POLYMEC_POLYGON_H

#include "core/point2.h"

/// \addtogroup geometry geometry
///@{

/// \class polygon
/// This class represents a polygon in the plane. Objects of this type are 
/// garbage-collected.
typedef struct polygon_t polygon_t;

/// Creates a new polygon in the plane given an ordered sequence of points 
/// that identify consecutive vertices.
polygon_t* polygon_new(point2_t* vertices, size_t num_vertices);

/// Creates a new polygon in the plane given a set of vertices and an 
/// ordering in which they are traversed.
polygon_t* polygon_new_with_ordering(point2_t* points, int* ordering, size_t num_points);

/// Creates a new convex polygon in the plane by applying the "gift-wrap" 
/// convex hull algorithm to the set of points.
polygon_t* polygon_giftwrap(point2_t* points, size_t num_points);

/// Creates a new nonconvex "star-shaped" polygon in the plane by ordering 
/// the given set of points according to their angles of displacement about 
/// a given center x0. The points must all be distinct from x0.
polygon_t* polygon_star(point2_t* x0, point2_t* points, size_t num_points);

/// Returns the number of vertices in the polygon.
size_t polygon_num_vertices(polygon_t* poly);

/// Returns the number of edges in the polygon (same as the number of 
/// vertices).
static inline size_t polygon_num_edges(polygon_t* poly)
{
  return polygon_num_vertices(poly);
}

/// Returns the ordering of the vertices in the polygon in terms of 
/// the original array of vertices given in the constructor.
int* polygon_ordering(polygon_t* poly);

/// Allows the traversal of the vertices in the polygon.
bool polygon_next_vertex(polygon_t* poly, int* pos, point2_t** vertex);

/// Returns the area of the polygon.
real_t polygon_area(polygon_t* poly);

/// Clones the polygon, returning an exact copy.
polygon_t* polygon_clone(polygon_t* poly);

/// Clips the given polygon by intersecting it with another. Here, poly is 
/// the polygon that is to be clipped, and it is modified in place.
/// This algorithm only works reliably if both polygons are convex.
void polygon_clip(polygon_t* poly, polygon_t* other);

#endif

