// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/tests/create_cubed_sphere.h"

//------------------------------------------------------------------------
// Coordinate mappings:
//   1. equiangular -> lat/lon (eq_to_ll)
//   2. lat/lon -> equiangular (
// These mappings are each other's inverse.
//------------------------------------------------------------------------
typedef struct
{
  int block;
  real_t R1, R2;
} equiangular_t;

static void equatorial_eq_to_ll_map_point(void* context, point_t* x, point_t* y)
{
  equiangular_t* eq = context;
  y->x = x->x + 0.5*M_PI*eq->block;
  y->y = atan(tan(x->y) * cos(x->x));
  y->z = eq->R1 + (eq->R2 - eq->R1) * x->z;
}

static void equatorial_eq_to_ll_J(void* context, point_t* x, real_t J[3][3])
{
  equiangular_t* eq = context;

  // Compute Gnomonic coordinates from our equiangular ones.
  real_t X = tan(x->x), Y = tan(x->y), delta2 = 1.0 + X*X + Y*Y;

  // Compute the change-of-basis matrix.
  J[0][0] = 1.0;
  J[0][1] = J[0][2] = 0.0;
  J[1][0] = -X*Y*sqrt(1.0+X*X)/delta2;
  J[1][1] = (1.0+Y*Y)*sqrt(1.0+X*X)/delta2;
  J[1][2] = 0.0;
  J[2][0] = J[2][1] = 0.0;
  J[2][2] = 1.0/(eq->R2 - eq->R1);
}

static void polar_eq_to_ll_map_point(void* context, point_t* x, point_t* y)
{
  equiangular_t* eq = context;

  // North or South?
  real_t s = (eq->block == 4) ? 1.0 : -1.0;

  real_t tan_x = tan(x->x), tan_y = tan(x->y);
  y->x = -atan2(tan_x, tan_y);
  y->y = s * atan(1.0/(sqrt(tan_x*tan_x + tan_y*tan_y)));
  y->z = eq->R1 + (eq->R2 - eq->R1) * x->z;
}

static void polar_eq_to_ll_J(void* context, point_t* x, real_t J[3][3])
{
  equiangular_t* eq = context;

  // North or South?
  real_t s = (eq->block == 4) ? 1.0 : -1.0;

  // Compute Gnomonic coordinates from our equiangular ones.
  real_t X = tan(x->x), Y = tan(x->y), delta2 = 1.0 + X*X + Y*Y;

  // Compute the change-of-basis matrix.
  J[0][0] = -s*Y*(1.0+X*X) / (X*X + Y*Y);
  J[0][1] = -s*X*(1.0+Y*Y) / (X*X + Y*Y);
  J[0][2] = 0.0;
  J[1][0] = -s*X*(1.0+X*X) / (delta2*sqrt(X*X + Y*Y));
  J[1][1] = -s*Y*(1.0+Y*Y) / (delta2*sqrt(X*X + Y*Y));
  J[1][2] = J[2][0] = J[2][1] = 0.0;
  J[2][2] = 1.0/(eq->R2 - eq->R1);
}

static void equatorial_ll_to_eq_map_point(void* context, point_t* x, point_t* y)
{
  equiangular_t* eq = context;
  y->x = x->x - 0.5*M_PI*eq->block;
  y->y = atan2(tan(x->y), cos(x->x));
  y->z = (x->z - eq->R1) / (eq->R2 - eq->R1);
}

static void equatorial_ll_to_eq_J(void* context, point_t* x, real_t J[3][3])
{
  equiangular_t* eq = context;

  // Compute Gnomonic coordinates from our equiangular ones.
  real_t X = tan(x->x), Y = tan(x->y), delta2 = 1.0 + X*X + Y*Y;

  // Compute the change-of-basis matrix.
  J[0][0] = 1.0;
  J[0][1] = J[0][2] = 0.0;
  J[1][0] = X*Y/(1.0 + Y*Y);
  J[1][1] = delta2 / ((1.0 + Y*Y) * sqrt(1.0 + X*X));
  J[1][2] = 0.0;
  J[2][0] = J[2][1] = 0.0;
  J[2][2] = eq->R2 - eq->R1;
}

static void polar_ll_to_eq_map_point(void* context, point_t* x, point_t* y)
{
  equiangular_t* eq = context;

  // North or South?
  real_t s = (eq->block == 4) ? 1.0 : -1.0;

  real_t tan_y = tan(x->y);
  y->x = -s * atan2(sin(x->x), tan_y);
  y->y = -atan2(cos(x->x), tan_y);
  y->z = (x->z - eq->R1) / (eq->R2 - eq->R1);
}

static void polar_ll_to_eq_J(void* context, point_t* x, real_t J[3][3])
{
  equiangular_t* eq = context;

  // North or South?
  real_t s = (eq->block == 4) ? 1.0 : -1.0;

  // Compute Gnomonic coordinates from our equiangular ones.
  real_t X = tan(x->x), Y = tan(x->y), delta2 = 1.0 + X*X + Y*Y;

  // Compute the change-of-basis matrix.
  J[0][0] = -s*Y / (1.0 + X*X);
  J[0][1] = -s*delta2*X / ((1.0+X*X)*sqrt(X*X + Y*Y));
  J[0][2] = 0.0;
  J[1][0] = s*X / (1.0 + Y*Y);
  J[1][1] = -s*delta2*Y / ((1.0+Y*Y)*sqrt(X*X + Y*Y));
  J[1][2] = 0.0;
  J[2][0] = J[2][1] = 0.0;
  J[2][2] = eq->R2 - eq->R1;
}

static coord_mapping_t* ll_to_eq(coord_mapping_t* eq_to_ll)
{
  equiangular_t* eq = coord_mapping_context(eq_to_ll);
  coord_mapping_vtable vtable = {.dtor = NULL};
  if (eq->block < 4)
  {
    vtable.map_point = equatorial_ll_to_eq_map_point;
    vtable.jacobian = equatorial_ll_to_eq_J;
  }
  else
  {
    vtable.map_point = polar_ll_to_eq_map_point;
    vtable.jacobian = polar_ll_to_eq_J;
  }

  char block_name[129];
  snprintf(block_name, 128, "lat/lon -> equiangular (block %d)", eq->block);
  coord_mapping_t* X = coord_mapping_new(block_name, eq, vtable);
  coord_mapping_set_inverse(X, eq_to_ll);
  retain_ref(eq_to_ll);
  return X;
}

//------------------------------------------------------------------------
//                        End coordinate mappings
//------------------------------------------------------------------------

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
                                 real_t R1, real_t R2)
{
  blockmesh_t* mesh = blockmesh_new(comm, patch_nxy, patch_nxy, patch_nz);

  // Add the four equatorial blocks.
  for (int b = 0; b < 4; ++b)
    blockmesh_add_block(mesh, block_nxy, block_nxy, block_nz);

  // Add the polar blocks.
  blockmesh_add_block(mesh, block_nxy, block_nxy, block_nz);
  blockmesh_add_block(mesh, block_nxy, block_nxy, block_nz);

  // Connect the east-west boundaries of the equatorial blocks.
  int east[4] = {1, 5, 6, 2};
  int west[4] = {0, 4, 7, 3};
  for (int b = 0; b < 4; ++b)
    blockmesh_connect_blocks(mesh, b, east, (b + 1) % 4, west);

  // Connect the equatorial blocks to the north pole block (4).
  int north[4][4] = {{2, 6, 7, 3},
                     {6, 7, 3, 2},
                     {7, 3, 2, 6},
                     {3, 2, 6, 7}};
  int north_block_boundaries[4][4] = {{1, 5, 4, 0},  // connection to block 0
                                      {2, 6, 5, 1},  // connection to block 1
                                      {3, 7, 6, 2},  // connection to block 2
                                      {0, 4, 7, 3}}; // connection to block 3
  for (int b = 0; b < 4; ++b)
    blockmesh_connect_blocks(mesh, b, north[b], 4, north_block_boundaries[b]);

  // Connect the equatorial blocks to the south pole block (5).
  int south[4][4] = {{1, 5, 4, 0},
                     {0, 1, 5, 4},
                     {4, 0, 1, 5},
                     {5, 4, 0, 1}};
  int south_block_boundaries[4][4] = {{2, 6, 7, 3},  // connection to block 0
                                      {1, 5, 6, 2},  // connection to block 1
                                      {0, 4, 5, 1},  // connection to block 2
                                      {3, 7, 4, 0}}; // connection to block 3
  for (int b = 0; b < 4; ++b)
    blockmesh_connect_blocks(mesh, b, south[b], 5, south_block_boundaries[b]);

  // Finalize the mesh and return it.
  blockmesh_finalize(mesh);
  return mesh;
}

coord_mapping_t* cubed_sphere_equator_block_coords(int block_index,
                                                   real_t R1,
                                                   real_t R2)
{
  ASSERT((block_index >= 0) && (block_index <= 3));
  ASSERT(R1 < R2);
  equiangular_t* eq = polymec_malloc(sizeof(equiangular_t));
  eq->block = block_index;
  eq->R1 = R1;
  eq->R2 = R2;
  coord_mapping_vtable vtable = {.map_point = equatorial_eq_to_ll_map_point,
                                 .jacobian = equatorial_eq_to_ll_J,
                                 .dtor = polymec_free};
  char block_name[129];
  snprintf(block_name, 128, "equatorial block %d", block_index);
  coord_mapping_t* X = coord_mapping_new(block_name, eq, vtable);

  coord_mapping_t* Xinv = ll_to_eq(X);
  coord_mapping_set_inverse(X, Xinv);
  return X;
}

coord_mapping_t* cubed_sphere_north_block_coords(real_t R1, real_t R2)
{
  ASSERT(R1 < R2);
  equiangular_t* eq = polymec_malloc(sizeof(equiangular_t));
  eq->block = 4;
  eq->R1 = R1;
  eq->R2 = R2;
  coord_mapping_vtable vtable = {.map_point = polar_eq_to_ll_map_point,
                                 .jacobian = polar_eq_to_ll_J,
                                 .dtor = polymec_free};
  coord_mapping_t* X = coord_mapping_new("north block", eq, vtable);

  coord_mapping_t* Xinv = ll_to_eq(X);
  coord_mapping_set_inverse(X, Xinv);
  return X;
}

coord_mapping_t* cubed_sphere_south_block_coords(real_t R1, real_t R2)
{
  ASSERT(R1 < R2);
  equiangular_t* eq = polymec_malloc(sizeof(equiangular_t));
  eq->block = 5;
  eq->R1 = R1;
  eq->R2 = R2;
  coord_mapping_vtable vtable = {.map_point = polar_eq_to_ll_map_point,
                                 .jacobian = polar_eq_to_ll_J,
                                 .dtor = polymec_free};
  coord_mapping_t* X = coord_mapping_new("south block", eq, vtable);

  coord_mapping_t* Xinv = ll_to_eq(X);
  coord_mapping_set_inverse(X, Xinv);
  return X;
}

