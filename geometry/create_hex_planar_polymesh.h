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

/// This function creates and returns a planar polymesh consisting of a set 
/// of hexagonal cells extending outward from a center cell.
/// \param [in] radius The number of hexes extending outward from the center hex. 
/// \param [in] h The length of the side of each hexagon.
/// \param [in] rotation The angle of rotation for the mesh (in radians), 
///                      relative to a reference hex with a flat bottom and top
///                      in the xy plane.
planar_polymesh_t* create_hex_planar_polymesh(size_t radius, 
                                              real_t h,
                                              real_t rotation);

///@}

#endif

