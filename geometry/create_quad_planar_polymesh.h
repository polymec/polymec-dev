// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_QUAD_PLANAR_POLYMESH_H
#define POLYMEC_CREATE_QUAD_PLANAR_POLYMESH_H

#include "core/point.h"
#include "geometry/planar_polymesh.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a planar polymesh consisting of nx x ny
/// quadrilateral cells spanning the given bounding box. The box can be
/// periodic in x, y, or both.
planar_polymesh_t* create_quad_planar_polymesh(size_t nx, size_t ny,
                                               bbox_t* bbox,
                                               bool periodic_in_x,
                                               bool periodic_in_y);

///@}

#endif

