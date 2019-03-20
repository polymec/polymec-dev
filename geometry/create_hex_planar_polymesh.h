// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
planar_polymesh_t* create_hex_planar_polymesh(int radius, real_t h);

///@}

#endif

