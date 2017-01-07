// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYHEDRON_H
#define POLYMEC_POLYHEDRON_H

#include "geometry/polygon.h"

// This class represents a polyhedron. Objects of this type are 
// garbage-collected.
typedef struct polyhedron_t polyhedron_t;

// Creates a new polyhedron given a sequence of (polygonal) faces 
// that consist of indices of vertices. vertices is an array of 
// length num_vertices that contains all of the vertices for the 
// polyhedron. faces is an array of vertex index pointers, of length
// num_faces, and faces[i][j] is the index of the jth vertex for face i.
// So the coordinates of the jth vertex of face i are stored in 
// vertices[faces[i][j]]. Accordingly, faces[i] contains 
// num_face_vertices[i] vertex indices.
polyhedron_t* polyhedron_new(point_t* vertices, int num_vertices,
                             int** faces, int* num_face_vertices, 
                             int num_faces);

// Returns the number of vertices in the polyhedron.
int polyhedron_num_vertices(polyhedron_t* poly);

// Returns the number of faces in the polyhedron.
int polyhedron_num_faces(polyhedron_t* poly);

// Allows the traversal of the vertices in the polyhedron.
bool polyhedron_next_vertex(polyhedron_t* poly, int* pos, point_t* vertex);

// Allows the traversal of the faces in the polyhedron.
bool polyhedron_next_face(polyhedron_t* poly, int* pos, polygon_t** face);

// Returns the volume of the polyhedron.
real_t polyhedron_volume(polyhedron_t* poly);

// Computes the centroid of the polyhedron, storing it in centroid.
void polyhedron_compute_centroid(polyhedron_t* poly, point_t* centroid);

// Clones the polyhedron, returning an exact copy.
polyhedron_t* polyhedron_clone(polyhedron_t* poly);

#endif

