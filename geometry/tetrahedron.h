// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TETRAHEDRON_H
#define POLYMEC_TETRAHEDRON_H

#include "core/point.h"

/// \addtogroup geometry geometry
///@{

/// \class tetrahedron
/// This class represents a tetrahedron.
/// \refcounted
typedef struct tetrahedron_t tetrahedron_t;

/// Creates a new tetrahedron whose vertices are set to a default
/// "reference state."
/// \memberof tetrahedron
tetrahedron_t* tetrahedron_new(void);

/// Sets the vertices (v1, v2, v3, v4) of the tetrahedron, recomputing its
/// geometric properties.
/// \memberof tetrahedron
void tetrahedron_set_vertices(tetrahedron_t* t,
                              point_t* v1,
                              point_t* v2,
                              point_t* v3,
                              point_t* v4);

/// Returns the volume of the tetrahedron.
/// \memberof tetrahedron
real_t tetrahedron_volume(tetrahedron_t* t);

/// Computes the centroid of the tetrahedron.
/// \memberof tetrahedron
void tetrahedron_compute_centroid(tetrahedron_t* t, point_t* centroid);

/// Computes the circumcenter of the tetrahedron. Note that this point may
/// lie outside the tetrahedron.
/// \memberof tetrahedron
void tetrahedron_compute_circumcenter(tetrahedron_t* t, point_t* circumcenter);

/// Returns true if the given point x falls within the tetrahedron.
/// \memberof tetrahedron
bool tetrahedron_contains_point(tetrahedron_t* t, point_t* x);

/// Returns true if the given point x falls inside the circumsphere of the tetrahedron.
/// \memberof tetrahedron
bool tetrahedron_circumsphere_contains_point(tetrahedron_t* t, point_t* x);

/// Returns true if the given point x falls on the surface of the circumsphere of the tetrahedron.
/// \memberof tetrahedron
bool tetrahedron_circumsphere_intersects_point(tetrahedron_t* t, point_t* x);

/// Computes the point y within or on the surface of the tetrahedron that is
/// closest to the given point x. If the point x is contained within the
/// tetrahedron, y is set to x.
/// \memberof tetrahedron
void tetrahedron_compute_nearest_point(tetrahedron_t* t, point_t* x, point_t* y);

/// Traverses the faces of the tetrahedron, returning true if a face remains
/// in the traversal and false if not. A face is retrieved in the sense that the
/// positions of its vertices are stored in v1, v2, and v3. Set *pos to 0 to
/// reset the traversal.
/// \memberof tetrahedron
bool tetrahedron_next_face(tetrahedron_t* t, int* pos, point_t* v1, point_t* v2, point_t* v3);

///@}

#endif

