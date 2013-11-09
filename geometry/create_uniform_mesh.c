// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/unordered_set.h"
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
  double Lx = bbox->x2 - bbox->x1;
  double Ly = bbox->y2 - bbox->y1;
  double Lz = bbox->z2 - bbox->z1;
  double dx = Lx/nx, dy = Ly/ny, dz = Lz/nz;

  // Create a uniform rectilinear mesh!
  double* xs = malloc(sizeof(double) * (nx+1));
  double* ys = malloc(sizeof(double) * (ny+1));
  double* zs = malloc(sizeof(double) * (nz+1));
  for (int i = 0; i <= nx; ++i)
    xs[i] = i*dx;
  for (int i = 0; i <= ny; ++i)
    ys[i] = i*dy;
  for (int i = 0; i <= nz; ++i)
    zs[i] = i*dz;

  mesh_t* mesh = create_rectilinear_mesh(comm, xs, nx+1, ys, ny+1, zs, nz+1);

  // Clean up.
  free(zs);
  free(ys);
  free(xs);

  return mesh;
}

