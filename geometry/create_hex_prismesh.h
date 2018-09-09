// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_HEX_PRISMESH_H
#define POLYMEC_CREATE_HEX_PRISMESH_H

#include "geometry/prismesh.h"
#include "geometry/hex_lattice.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a prismesh consisting of a set of 
/// hexagonal prism cells with the given number of cells extending outward 
/// from a center hex, and with the given alignment. 
/// \param comm [in] The MPI communicator on which the distributed prismesh 
///                  is defined.
/// \param alignment [in] The alignment of the hexes in the mesh (x- or y-).
/// \param radius [in] The number of hexes extending outward from the center hex. 
/// \param h [in] The length of the side of each hexagon.
/// \param nz [in] The number of cells in the mesh along the z axis.
/// \param z1 [in] The z coordinate of the bottom of the mesh.
/// \param z2 [in] The z coordinate of the top of the mesh.
/// \param periodic_in_z [in] True if the mesh is periodic along the z axis, 
///                           false if not.
prismesh_t* create_hex_prismesh(MPI_Comm comm,
                                hex_lattice_align_t alignment,
                                size_t radius, real_t h,
                                size_t nz, real_t z1, real_t z2,
                                bool periodic_in_z);

///@}

#endif

