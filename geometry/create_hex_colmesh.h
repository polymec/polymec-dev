// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_HEX_COLMESH_H
#define POLYMEC_CREATE_HEX_COLMESH_H

#include "geometry/colmesh.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a colmesh consisting of a set of 
/// hexagonal prism cells with the given number of cells extending outward 
/// from a center hex, and with the given alignment. 
/// \param [in] comm The MPI communicator on which the distributed colmesh 
///                  is defined.
/// \param [in] radius The number of hexes extending outward from the center hex. 
/// \param [in] h  The length of the side of each hexagon.
/// \param [in] nz The number of cells in the mesh along the z axis.
/// \param [in] z1 The z coordinate of the bottom of the mesh.
/// \param [in] z2 The z coordinate of the top of the mesh.
/// \param [in] periodic_in_z True if the mesh is periodic along the z axis, 
///                           false if not.
colmesh_t* create_hex_colmesh(MPI_Comm comm, 
                              size_t radius, real_t h, 
                              size_t nz, real_t z1, real_t z2,
                              bool periodic_in_z);

///@}

#endif

