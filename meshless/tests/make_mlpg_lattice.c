// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_point_cloud.h"
#include "core/st_func.h"
#include "geometry/create_point_lattice.h"
#include "make_mlpg_lattice.h"

static void n_eval(void* context, point_t* x, real_t* n)
{
  bbox_t* bbox = context;
  n[0] = n[1] = n[2] = 0.0;
  if (x->x <= bbox->x1)
    n[0] = -1.0;
  else if (x->x >= bbox->x2)
    n[0] =  1.0;
  else if (x->y <= bbox->y1)
    n[1] = -1.0;
  else if (x->y >= bbox->y2)
    n[1] =  1.0;
  else if (x->z <= bbox->z1)
    n[2] = -1.0;
  else if (x->z >= bbox->z2)
    n[2] =  1.0;
}

void make_mlpg_lattice(bbox_t* bbox, int nx, int ny, int nz, real_t R_over_dx,
                       point_cloud_t** domain,
                       real_t** extents,
                       stencil_t** neighborhoods)
{
  real_t dx = (bbox->x2 - bbox->x1)/nx, 
         dy = (bbox->y2 - bbox->y1)/ny, 
         dz = (bbox->z2 - bbox->z1)/nz;
  if (nx == 1)
    dx = MIN(dy, dz);
  else if (ny == 1)
    dy = MIN(dx, dz);
  else if (nz == 1)
    dz = MIN(dx, dy);

  // Set up a lattice for interior points.
  *domain = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, nz, bbox);

  // Add boundary points for all the faces of the domain.
  point_cloud_t* bcloud;
  bbox->x1 = bbox->x2 = 0.0;
  bcloud = create_uniform_point_lattice(MPI_COMM_SELF, 1, ny, nz, bbox);
  point_cloud_unite(*domain, bcloud, "-x");
  point_cloud_free(bcloud);
  bbox->x1 = bbox->x2 = nx*dx;
  bcloud = create_uniform_point_lattice(MPI_COMM_SELF, 1, ny, nz, bbox);
  point_cloud_unite(*domain, bcloud, "+x");
  point_cloud_free(bcloud);
  bbox->x1 = 0.0, bbox->x2 = nx*dx;
  bbox->y1 = bbox->y2 = 0.0;
  bcloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, 1, nz, bbox);
  point_cloud_unite(*domain, bcloud, "-y");
  point_cloud_free(bcloud);
  bbox->y1 = bbox->y2 = ny*dy;
  bcloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, 1, nz, bbox);
  point_cloud_unite(*domain, bcloud, "+y");
  point_cloud_free(bcloud);
  bbox->y1 = 0.0, bbox->y2 = ny*dy;
  bbox->z1 = bbox->z2 = 0.0;
  bcloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, 1, bbox);
  point_cloud_unite(*domain, bcloud, "-z");
  point_cloud_free(bcloud);
  bbox->z1 = bbox->z2 = nz*dz;
  bcloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, 1, bbox);
  point_cloud_unite(*domain, bcloud, "+z");
  point_cloud_free(bcloud);

  // Make an omni-boundary tag.
  tagger_unite_tag((*domain)->tags, "boundary", "-x");
  tagger_unite_tag((*domain)->tags, "boundary", "+x");
  tagger_unite_tag((*domain)->tags, "boundary", "-y");
  tagger_unite_tag((*domain)->tags, "boundary", "+y");
  tagger_unite_tag((*domain)->tags, "boundary", "-z");
  tagger_unite_tag((*domain)->tags, "boundary", "+z");

  // Do partitioning.
  exchanger_t* ex = partition_point_cloud(domain, MPI_COMM_WORLD, NULL, 1.05);
  exchanger_free(ex);

  // Set up a normal vector field for the boundary of the domain.
  point_t x0 = {.x = 0.5, .y = 0.5, .z = 0.5/nx};
  sp_func_vtable n_vtable = {.eval = n_eval};
  sp_func_t* n = sp_func_new("mlpg lattice normal", bbox, n_vtable,
                             SP_FUNC_INHOMOGENEOUS, 3);
  point_cloud_set_property(*domain, "normal", n, NULL);

  // Set up a "radius" field to measure point extents.
  int num_local_points = (*domain)->num_points;
  int num_remote_points = (*domain)->num_ghosts;
  int num_points = num_local_points + num_remote_points;
  *extents = polymec_malloc(sizeof(real_t) * num_points);
  for (int i = 0; i < num_local_points; ++i)
    (*extents)[i] = R_over_dx * dx;

  // Create the stencil if requested.
  if (neighborhoods != NULL)
    *neighborhoods = distance_based_point_stencil_new(*domain, *extents);
}


