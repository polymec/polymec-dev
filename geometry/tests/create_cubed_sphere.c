// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/blockmesh.h"

typedef struct
{
  int block;
  real_t R1, R2;
} equiangular_t;

static void equatorial_map_point(void* context, point_t* x, point_t* y)
{
  equiangular_t* eq = context;
  y->x = x->x + 0.5*M_PI*eq->block;
  y->y = atan(tan(x->y) * cos(x->x));
  y->z = eq->R1 + (eq->R2 - eq->R1) * x->z;
}

static void equatorial_J(void* context, point_t* x, tensor2_t* J)
{
  equiangular_t* eq = context;

  // Compute Gnomonic coordinates from our equiangular ones.
  real_t X = tan(x->x), Y = tan(x->y), delta2 = 1.0 + X*X + Y*Y;

  // Compute the change-of-basis matrix.
  J->xx = 1.0; 
  J->xy = J->xz = 0.0;
  J->yx = X*Y/(1.0 + Y*Y);
  J->yy = delta2 / ((1.0 + Y*Y) * sqrt(1.0 + X*X));
  J->yz = 0.0;
  J->zx = J->zy = 0.0;
  J->zz = eq->R2 - eq->R1;
}

static void polar_map_point(void* context, point_t* x, point_t* y)
{
  equiangular_t* eq = context;

  // North or South?
  real_t s = (eq->block == 4) ? 1.0 : -1.0;

  y->x = -atan2(tan(x->x), tan(x->y));
  real_t tan_x = tan(x->x), tan_y = tan(x->y);
  y->y = s * atan(1.0/(sqrt(tan_x*tan_x + tan_y*tan_y)));
  y->z = eq->R1 + (eq->R2 - eq->R1) * x->z;
}

static void polar_J(void* context, point_t* x, tensor2_t* J)
{
  equiangular_t* eq = context;

  // North or South?
  real_t s = (eq->block == 4) ? 1.0 : -1.0;

  // Compute Gnomonic coordinates from our equiangular ones.
  real_t X = tan(x->x), Y = tan(x->y), delta2 = 1.0 + X*X + Y*Y;

  // Compute the change-of-basis matrix.
  J->xx = -s*Y / (1.0 + X*X);
  J->xy = -s*delta2*X / ((1.0+X*X)*sqrt(X*X + Y*Y));
  J->yx = s*X / (1.0 + Y*Y);
  J->yy = -s*delta2*Y / ((1.0+Y*Y)*sqrt(X*X + Y*Y));
  J->yz = 0.0;
  J->zx = J->zy = 0.0;
  J->zz = eq->R2 - eq->R1;
}

static coord_mapping_t* create_equator_block_coords(int block_index,
                                                    real_t R1,
                                                    real_t R2)
{
  ASSERT((block_index >= 0) && (block_index <= 3));
  ASSERT(R1 < R2);
  equiangular_t* eq = polymec_malloc(sizeof(equiangular_t));
  eq->block = block_index;
  eq->R1 = R1;
  eq->R2 = R2;
  coord_mapping_vtable vtable = {.map_point = equatorial_map_point,
                                 .jacobian = equatorial_J,
                                 .dtor = polymec_free};
  char block_name[129];
  snprintf(block_name, 128, "equatorial block %d", block_index);
  return coord_mapping_new(block_name, eq, vtable);
}

static coord_mapping_t* create_north_block_coords(real_t R1, real_t R2)
{
  ASSERT(R1 < R2);
  equiangular_t* eq = polymec_malloc(sizeof(equiangular_t));
  eq->block = 4;
  eq->R1 = R1;
  eq->R2 = R2;
  coord_mapping_vtable vtable = {.map_point = polar_map_point,
                                 .jacobian = polar_J,
                                 .dtor = polymec_free};
  return coord_mapping_new("north block", eq, vtable);
}

static coord_mapping_t* create_south_block_coords(real_t R1, real_t R2)
{
  ASSERT(R1 < R2);
  equiangular_t* eq = polymec_malloc(sizeof(equiangular_t));
  eq->block = 5;
  eq->R1 = R1;
  eq->R2 = R2;
  coord_mapping_vtable vtable = {.map_point = polar_map_point,
                                 .jacobian = polar_J,
                                 .dtor = polymec_free};
  return coord_mapping_new("south block", eq, vtable);
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

  // We use equiangular coordinates, which span [-pi/4, pi/4] in the
  // first two coordinates.
  bbox_t domain = {.x1 = -0.25 * M_PI, .x2 = 0.25 * M_PI,
                   .y1 = -0.25 * M_PI, .y2 = 0.25 * M_PI,
                   .z1 = 0.0, .z2 = 1.0};

  // Add the four equatorial blocks.
  for (int b = 0; b < 4; ++b)
  {
    blockmesh_add_block(mesh, &domain, create_equator_block_coords(b, R1, R2),
                        block_nxy, block_nxy, block_nz);
  }

  // Add the polar blocks.
  blockmesh_add_block(mesh, &domain, create_north_block_coords(R1, R2),
                      block_nxy, block_nxy, block_nz);
  blockmesh_add_block(mesh, &domain, create_south_block_coords(R1, R2),
                      block_nxy, block_nxy, block_nz);

  // Connect the east-west boundaries of the equatorial blocks.
  int east[4] = {1, 5, 6, 2};
  int west[4] = {0, 4, 7, 3};
  for (int b = 0; b < 4; ++b)
    blockmesh_connect_blocks(mesh, b, east, (b + 1) % 4, west);

  // Connect the equatorial blocks to the north pole block (4).
  int north[4] = {2, 6, 7, 3};
  int north_block_boundaries[4][4] = {{1, 5, 4, 0},  // connection to block 0
                                      {2, 6, 5, 1},  // connection to block 1
                                      {3, 7, 6, 2},  // connection to block 2
                                      {0, 4, 7, 3}}; // connection to block 3
  for (int b = 0; b < 4; ++b)
    blockmesh_connect_blocks(mesh, b, north, 4, north_block_boundaries[b]);

  // Connect the equatorial blocks to the south pole block (5).
  int south[4] = {1, 5, 4, 0};
  int south_block_boundaries[4][4] = {{2, 6, 7, 3},  // connection to block 0
                                      {1, 5, 6, 2},  // connection to block 1
                                      {0, 4, 5, 1},  // connection to block 2
                                      {3, 7, 4, 0}}; // connection to block 3
  for (int b = 0; b < 4; ++b)
    blockmesh_connect_blocks(mesh, b, south, 5, south_block_boundaries[b]);

  // Finalize the mesh and return it.
  blockmesh_finalize(mesh);
  return mesh;
}

