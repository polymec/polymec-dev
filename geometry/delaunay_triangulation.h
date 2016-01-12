// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

