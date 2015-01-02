// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_DELAUNAY_TRIANGULATION_H
#define POLYMEC_DELAUNAY_TRIANGULATION_H

#include "core/point.h"

// This class represents a Delaunay triangulation in 3D.
typedef struct delaunay_triangulation_t delaunay_triangulation_t;

// Creates a new Delaunay triangulation from the given set of points.
delaunay_triangulation_t* delaunay_triangulation_new(point_t* points, int num_points);

// Frees the given triangulation.
void delaunay_triangulation_free(delaunay_triangulation_t* t);

// Returns the number of vertices in the triangulation.
int delaunay_triangulation_num_vertices(delaunay_triangulation_t* t);

// Retrieves the coordinates of the vertices with the given indices from 
// the triangulation.
void delaunay_triangulation_get_vertices(delaunay_triangulation_t* t, int* indices, int num_vertices, point_t* vertices);

// Returns the number of tetrahedra in the triangulation.
int delaunay_triangulation_num_tetrahedra(delaunay_triangulation_t* t);

// Allows traversal over each tetrahedron in the triangulation, storing 
// the indices of the vertices in v1, v2, v3, v4. Set *pos to 0 
// to reset the iteration.
bool delaunay_triangulation_next(delaunay_triangulation_t* t,
                                 int* pos, int* v1, int* v2, int* v3, int* v4);

#endif

