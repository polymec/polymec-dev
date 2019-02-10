// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/blockmesh.h"

static coord_mapping_t* create_equator_block_coords(int block_index,
                                                    real_t R1,
                                                    real_t R2)
{
  return NULL;
}

static coord_mapping_t* create_north_block_coords(real_t R1, real_t R2)
{
  return NULL;
}

static coord_mapping_t* create_south_block_coords(real_t R1, real_t R2)
{
  return NULL;
}

// This function creates a cubed-sphere blockmesh for testing purposes.
// This mesh stitches together 6 blocks to form a spherical shell with 
// inner radius R1 and outer radius R2. 
// * Each block has block_nxy patches in the azimuthal and colatitudinal 
//   directions, and block_nz patches in the radial direction.
// * Each patch has patch_nxy cells in the azimuthal and colatitudinal
//   directions, and patch_nz cells in the radial direction.
blockmesh_t* create_cubed_sphere(MPI_Comm comm,
                                 int block_nxy, int block_nz,
                                 int patch_nxy, int patch_nz,
                                 real_t R1, real_t R2);
blockmesh_t* create_cubed_sphere(MPI_Comm comm,
                                 int block_nxy, int block_nz,
                                 int patch_nxy, int patch_nz,
                                 real_t R1, real_t R2)
{
  blockmesh_t* mesh = blockmesh_new(comm, patch_nxy, patch_nxy, patch_nz);

  // Add the four equatorial blocks.
  for (int b = 0; b < 4; ++b)
  {
    blockmesh_add_block(mesh, create_equator_block_coords(b, R1, R2),
                        block_nxy, block_nxy, block_nz);
  }

  // Add the polar blocks.
  blockmesh_add_block(mesh, create_north_block_coords(R1, R2),
                      block_nxy, block_nxy, block_nz);
  blockmesh_add_block(mesh, create_south_block_coords(R1, R2),
                      block_nxy, block_nxy, block_nz);

  // Connect the blocks.

  // Finalize the mesh and return it.
  blockmesh_finalize(mesh);
  return mesh;
}

