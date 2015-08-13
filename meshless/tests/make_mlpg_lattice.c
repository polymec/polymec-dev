// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_point_cloud.h"
#include "geometry/create_point_lattice.h"
#include "make_mlpg_lattice.h"

void make_mlpg_lattice(int nx, int ny, int nz, real_t R_over_dx,
                       point_cloud_t** domain,
                       real_t** extents,
                       stencil_t** neighborhoods)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  *domain = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, nz, &bbox);
  exchanger_t* ex = partition_point_cloud(domain, MPI_COMM_WORLD, NULL, 1.05);
  exchanger_free(ex);
  real_t dx = 1.0 / nx;

  // Set up a "radius" field to measure point extents.
  int num_local_points = (*domain)->num_points;
  int num_remote_points = (*domain)->num_ghosts;
  int num_points = num_local_points + num_remote_points;
  *extents = polymec_malloc(sizeof(real_t) * num_points);
  for (int i = 0; i < num_local_points; ++i)
    (*extents)[i] = R_over_dx * dx;

  // Tag boundary points.
  int_array_t* bpoints = int_array_new();
  for (int i = 0; i < num_local_points; ++i)
  {
    point_t* x = &((*domain)->points[i]);
    if (((nx > 1) && ((x->x == 0.0) || (x->x == 1.0))) || 
        ((ny > 1) && ((x->y == 0.0) || (x->y == 1.0))) || 
        ((nz > 1) && ((x->z == 0.0) || (x->z == 1.0))))
      int_array_append(bpoints, i);
  }
  int* btag = point_cloud_create_tag(*domain, "boundary", bpoints->size);
  memcpy(btag, bpoints->data, sizeof(int) * bpoints->size);
  int_array_free(bpoints);

  // Create the stencil if requested.
  if (neighborhoods != NULL)
    *neighborhoods = distance_based_point_stencil_new(*domain, *extents);
}


