// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "geometry/create_point_lattice.h"
#include "geometry/cubic_lattice.h"

point_cloud_t* create_point_lattice(MPI_Comm comm, 
                                    real_t* xs, int nxs, 
                                    real_t* ys, int nys, 
                                    real_t* zs, int nzs)
{
  ASSERT(xs != NULL);
  ASSERT(nxs > 0);
  ASSERT(ys != NULL);
  ASSERT(nys > 0);
  ASSERT(zs != NULL);
  ASSERT(nzs > 0);

  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Create a new point cloud.
  int points_per_proc = nxs * nys * nzs / nprocs;
  int num_local_points = (rank < (nprocs-1)) ? points_per_proc
                                             : nxs * nys * nzs - rank * points_per_proc;
  point_cloud_t* cloud = point_cloud_new(comm, num_local_points);

  // Create a cubic lattice object for indexing.
  cubic_lattice_t* lattice = cubic_lattice_new(nxs, nys, nzs);

  // Now create the points in the cloud.
  for (int p = 0; p < num_local_points; ++p)
  {
    // Figure out (i, j, k) indices for this cell.
    uint64_t global_index = rank*points_per_proc + p, i, j, k;
    cubic_lattice_get_cell_triple(lattice, global_index, &i, &j, &k);

    // Set the point's coordinates.
    point_t* xp = &cloud->points[p];
    xp->x = xs[i];
    xp->y = ys[j];
    xp->z = zs[k];
  }

  return cloud;
}

point_cloud_t* create_uniform_point_lattice(MPI_Comm comm, int nx, int ny, int nz, bbox_t* bbox)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(bbox != NULL);

  ASSERT(bbox->x2 > bbox->x1);
  ASSERT(bbox->y2 > bbox->y1);
  ASSERT(bbox->z2 > bbox->z1);

  // Lattice spacings.
  real_t Lx = bbox->x2 - bbox->x1;
  real_t Ly = bbox->y2 - bbox->y1;
  real_t Lz = bbox->z2 - bbox->z1;
  real_t dx = (nx == 1) ? 0.0 : Lx/nx, 
         dy = (ny == 1) ? 0.0 : Ly/ny, 
         dz = (nz == 1) ? 0.0 : Lz/nz;

  // Create a uniform rectilinear lattice!
  real_t* xs = polymec_malloc(sizeof(real_t) * nx);
  real_t* ys = polymec_malloc(sizeof(real_t) * ny);
  real_t* zs = polymec_malloc(sizeof(real_t) * nz);
  for (int i = 0; i < nx; ++i)
    xs[i] = bbox->x1 + (0.5+i)*dx;
  for (int i = 0; i < ny; ++i)
    ys[i] = bbox->y1 + (0.5+i)*dy;
  for (int i = 0; i < nz; ++i)
    zs[i] = bbox->z1 + (0.5+i)*dz;

  point_cloud_t* cloud = create_point_lattice(comm, xs, nx, ys, ny, zs, nz);

  // Clean up.
  polymec_free(zs);
  polymec_free(ys);
  polymec_free(xs);

  return cloud;
}

