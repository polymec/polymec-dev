// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_uniform_polymesh.h"
#include "geometry/create_rectilinear_polymesh.h"
#include "geometry/cubic_lattice.h"

polymesh_t* create_uniform_polymesh(MPI_Comm comm, int nx, int ny, int nz, bbox_t* bbox)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(bbox != NULL);

  ASSERT(bbox->x2 > bbox->x1);
  ASSERT(bbox->y2 > bbox->y1);
  ASSERT(bbox->z2 > bbox->z1);

  // Grid spacings.
  real_t Lx = bbox->x2 - bbox->x1;
  real_t Ly = bbox->y2 - bbox->y1;
  real_t Lz = bbox->z2 - bbox->z1;
  real_t dx = Lx/nx, dy = Ly/ny, dz = Lz/nz;

  // Create a uniform rectilinear mesh!
  real_t* xs = polymec_malloc(sizeof(real_t) * (nx+1));
  real_t* ys = polymec_malloc(sizeof(real_t) * (ny+1));
  real_t* zs = polymec_malloc(sizeof(real_t) * (nz+1));
  for (int i = 0; i <= nx; ++i)
    xs[i] = bbox->x1 + i*dx;
  for (int i = 0; i <= ny; ++i)
    ys[i] = bbox->y1 + i*dy;
  for (int i = 0; i <= nz; ++i)
    zs[i] = bbox->z1 + i*dz;

  polymesh_t* mesh = create_rectilinear_polymesh(comm, xs, nx+1, ys, ny+1, zs, nz+1);

  // Clean up.
  polymec_free(zs);
  polymec_free(ys);
  polymec_free(xs);

  return mesh;
}

polymesh_t* create_uniform_polymesh_on_rank(MPI_Comm comm, int rank,
                                            int nx, int ny, int nz, bbox_t* bbox)
{
  ASSERT(comm != MPI_COMM_SELF);
  ASSERT(rank >= 0);

  int my_rank, nprocs;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_rank(comm, &nprocs);

  if (rank > nprocs)
    polymec_error("create_uniform_polymesh_on_rank: invalid rank: %d", rank);

  polymesh_t* mesh = NULL;
  if (my_rank == rank)
    mesh = create_uniform_polymesh(MPI_COMM_SELF, nx, ny, nz, bbox);
  else
  {
    // Initialize serializers.
    serializer_t* s = cubic_lattice_serializer();
    s = bbox_serializer();
  }
  return mesh;
}
