// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_HEX_PLANAR_POLYMESH_H
#define POLYMEC_CREATE_HEX_PLANAR_POLYMESH_H

#include "core/point.h"
#include "geometry/planar_polymesh.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a planar polymesh consisting of nx x ny 
/// regular hexagonal cells filling the given bounding box. The cells on the 
/// boundary are cropped to fit the bounding box unless the box is periodic
/// on that boundary.
planar_polymesh_t* create_hex_planar_polymesh(size_t nx, size_t ny, 
                                              bbox_t* bbox,
                                              bool periodic_in_x,
                                              bool periodic_in_y);

///@}

#endif

