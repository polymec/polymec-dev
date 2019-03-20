// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
/// This class represents a polygon in the plane.
/// \refcounted
typedef struct polygon_t polygon_t;

/// Creates a new polygon in the plane given an ordered sequence of points
/// that identify consecutive vertices.
/// \memberof polygon
polygon_t* polygon_new(point2_t* vertices, size_t num_vertices);

/// Creates a new polygon in the plane given a set of vertices and an
/// ordering in which they are traversed.
/// \memberof polygon
polygon_t* polygon_new_with_ordering(point2_t* vertices, int* ordering, size_t num_vertices);

/// Creates a new convex polygon in the plane by applying the "gift-wrap"
/// convex hull algorithm to the set of points.
/// \memberof polygon
polygon_t* polygon_giftwrap(point2_t* points, size_t num_points);

/// Creates a new nonconvex "star-shaped" polygon in the plane by ordering
/// the given set of points according to their angles of displacement about
/// a given center x0. The points must all be distinct from x0.
/// \memberof polygon
polygon_t* polygon_star(point2_t* x0, point2_t* points, size_t num_points);

/// Returns the number of vertices in the polygon.
/// \memberof polygon
size_t polygon_num_vertices(polygon_t* poly);

/// Sets the vertices for an existing polygon.
/// \param [in] vertices An array of 2D points representing the new vertices of the
///                      polygon.
/// \param [in] num_vertices The number of vertices.
/// \memberof polygon
void polygon_set_vertices(polygon_t* poly, point2_t* vertices, size_t num_vertices);

/// Returns the number of edges in the polygon (same as the number of
/// vertices).
/// \memberof polygon
static inline size_t polygon_num_edges(polygon_t* poly)
{
  return polygon_num_vertices(poly);
}

/// Allows the traversal of the vertices in the polygon.
/// \memberof polygon
bool polygon_next_vertex(polygon_t* poly, int* pos, point2_t* vertex);

/// Returns the area of the polygon.
/// \memberof polygon
real_t polygon_area(polygon_t* poly);

/// Computes the centroid of the polygon.
/// \memberof polygon
void polygon_compute_centroid(polygon_t* poly, point2_t* centroid);

/// Clones the polygon, returning an exact copy.
/// \memberof polygon
polygon_t* polygon_clone(polygon_t* poly);

///@}

#endif

