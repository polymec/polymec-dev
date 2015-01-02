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

#include "geometry/create_uniform_mesh.h"
#include "geometry/create_rectilinear_mesh.h"
#include "geometry/cubic_lattice.h"

mesh_t* create_uniform_mesh(MPI_Comm comm, int nx, int ny, int nz, bbox_t* bbox)
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

  mesh_t* mesh = create_rectilinear_mesh(comm, xs, nx+1, ys, ny+1, zs, nz+1);

  // Clean up.
  polymec_free(zs);
  polymec_free(ys);
  polymec_free(xs);

  return mesh;
}

