// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/blockmesh.h"
#include "geometry/coord_mapping.h"

// This function creates a ring of 4 equatorial blocks for testing purposes.
// * Each block has block_nxy patches in the azimuthal and colatitudinal
//   directions, and block_nz patches in the radial direction.
// * Each patch has patch_nxy cells in the azimuthal and colatitudinal
//   directions, and patch_nz cells in the radial direction.
blockmesh_t* create_multiblock_mesh(MPI_Comm comm,
                                    int block_nxy, int block_nz,
                                    int patch_nxy, int patch_nz,
                                    real_t R1, real_t R2);

// Here's a coordinate mapping for equatorial blocks (0-3) on a cubed sphere.
coord_mapping_t* cubed_sphere_equator_block_coords(int block_index,
                                                   real_t R1,
                                                   real_t R2);

