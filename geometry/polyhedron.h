// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYHEDRON_H
#define POLYMEC_POLYHEDRON_H

#include "core/point.h"

/// \addtogroup geometry geometry
///@{

/// \class polyhedron
/// A polyhedron in 3D space.
/// \refcounted
typedef struct polyhedron_t polyhedron_t;

/// Creates a new polyhedron given a sequence of (polygonal) faces. These
/// faces are represented by ordered lists of vertices.
/// \param [in] vertices An array of length num_vertices containing all of
///                      the vertices for the polyhedron.
/// \param [in] num_vertices The number of distinct vertices for all of the
///                          polygonal faces in the polyhedron.
/// \param [in] faces An array of vertex index pointers representing ordered
///                   lists of vertices. Specifically, faces[i] is the list of
///                   vertices for face i, traversed counterclockwise.
///                   faces[i][j] is the index of the jth vertex for face i.
///                   So the coordinates of the jth vertex of face i are
///                   stored in vertices[faces[i][j]].
/// \param num_face_vertices [in] An array of length num_faces giving the number
///                               of vertices in each face. num_face_vertices[i]
///                               gives the number of vertices in the ith face.
/// \param num_faces [in] The number of faces in the polyhedron.
/// \returns A newly created polyhedron.
/// \memberof polyhedron
polyhedron_t* polyhedron_new(point_t* vertices, size_t num_vertices,
                             int** faces, size_t* num_face_vertices,
                             size_t num_faces);

/// Returns the number of vertices in the polyhedron.
/// \memberof polyhedron
size_t polyhedron_num_vertices(polyhedron_t* poly);

/// Returns the number of faces in the polyhedron.
/// \memberof polyhedron
size_t polyhedron_num_faces(polyhedron_t* poly);

/// Allows the traversal of the vertices in the polyhedron. These vertices
/// are not traversed in any meaningful order.
/// \memberof polyhedron
bool polyhedron_next_vertex(polyhedron_t* poly, int* pos, point_t* vertex);

/// Allows the traversal of the faces in the polyhedron.
/// \param pos [in,out] Controls the traversal. Set to 0 to reset the traversal.
/// \param face_vertices [out] Stores an array of vertices for the next face.
/// \param num_face_vertices [out] Stores the number of vertices in the next face.
/// \param face_area [out] If non-NULL, stores the area of the next face.
/// \param face_normal [out] If non-NULL, stores the normal vector for the next face.
/// \memberof polyhedron
bool polyhedron_next_face(polyhedron_t* poly, int* pos,
                          point_t** face_vertices,
                          size_t* num_face_vertices,
                          real_t* face_area,
                          vector_t* face_normal);

/// Returns the volume of the polyhedron.
/// \memberof polyhedron
real_t polyhedron_volume(polyhedron_t* poly);

/// Computes the centroid of the polyhedron, storing it in centroid.
/// \memberof polyhedron
void polyhedron_compute_centroid(polyhedron_t* poly, point_t* centroid);

/// Clones the polyhedron, returning an exact copy.
/// \memberof polyhedron
polyhedron_t* polyhedron_clone(polyhedron_t* poly);

///@}

#endif

