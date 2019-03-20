// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PLANE_SD_FUNC_H
#define POLYMEC_PLANE_SD_FUNC_H

#include "geometry/sd_func.h"
#include "core/point2.h"

/// \addtogroup geometry geometry
///@{

/// This signed distance function represents a plane with a given normal
/// vector and point.
/// \relates sd_func
sd_func_t* plane_sd_func_new(vector_t* n, point_t* x);

/// Construct a plane that contains the three points p1, p2, and p3,
/// assuming that these points are not co-linear.
/// \relates sd_func
sd_func_t* plane_sd_func_from_points(point_t* p1, point_t* p2, point_t* p3);

/// Projects the given point x to a 2D point y within the plane.
/// \relates sd_func
void plane_sd_func_project(sd_func_t* plane, point_t* x, point2_t* y);

/// Maps a point y in the plane to a 3D point x. This is the inverse operation
/// for \ref plane_sd_func_project.
/// \relates sd_func
void plane_sd_func_embed(sd_func_t* plane, point2_t* y, point_t* x);

///@}

#endif

