// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
// a given center x0.
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

