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

/// \enum hex_lattice_align_t
/// This type identifies two possible alignments for hexagonal lattices.
typedef enum
{
  /// An x-aligned lattice has edges separating cells along the x axis.
  HEX_LATTICE_X_ALIGNED = 0,
  /// A y-aligned lattice has edges separating cells to the y axis.
  HEX_LATTICE_Y_ALIGNED = 1
} hex_lattice_align_t;

/// This function creates and returns a planar polymesh consisting of a set 
/// of hexagonal cells extending outward from a center cell.
/// \param alignment [in] The alignment of the hexes in the mesh (x- or y-).
/// \param radius [in] The number of hexes extending outward from the center hex. 
/// \param h [in] The length of the side of each hexagon.
planar_polymesh_t* create_hex_planar_polymesh(hex_lattice_align_t alignment,
                                              size_t radius, real_t h);

///@}

#endif

