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

#ifndef POLYMEC_TETRAHEDRON_H
#define POLYMEC_TETRAHEDRON_H

#include "core/point.h"

// This class represents a tetrahedron.
// Objects of this type are garbage-collected.
typedef struct tetrahedron_t tetrahedron_t;

// Creates a new tetrahedron whose vertices are set to a default 
// "reference state."
tetrahedron_t* tetrahedron_new();

// Sets the vertices (v1, v2, v3, v4) of the tetrahedron, recomputing its 
// geometric properties.
void tetrahedron_set_vertices(tetrahedron_t* t, 
                              point_t* v1,
                              point_t* v2,
                              point_t* v3,
                              point_t* v4);

// Returns the volume of the tetrahedron.
real_t tetrahedron_volume(tetrahedron_t* t);

// Computes the centroid of the tetrahedron.
void tetrahedron_compute_centroid(tetrahedron_t* t, point_t* centroid);

// Computes the circumcenter of the tetrahedron. Note that this point may 
// lie outside the tetrahedron.
void tetrahedron_compute_circumcenter(tetrahedron_t* t, point_t* circumcenter);

// Returns true if the given point x falls within the tetrahedron.
bool tetrahedron_contains_point(tetrahedron_t* t, point_t* x);

// Returns true if the given point x falls within the circumsphere of the tetrahedron.
bool tetrahedron_circumsphere_contains_point(tetrahedron_t* t, point_t* x);

// Computes the point y within or on the surface of the tetrahedron that is 
// closest to the given point x. If the point x is contained within the 
// tetrahedron, y is set to x.
void tetrahedron_compute_nearest_point(tetrahedron_t* t, point_t* x, point_t* y);

// Traverses the faces of the tetrahedron, returning true if a face remains
// in the traversal and false if not. A face is retrieved in the sense that the 
// positions of its vertices are stored in v1, v2, and v3. Set *pos to 0 to 
// reset the traversal.
bool tetrahedron_next_face(tetrahedron_t* t, int* pos, point_t* v1, point_t* v2, point_t* v3);

#endif

