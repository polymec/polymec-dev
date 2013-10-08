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

#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/create_voronoi_mesh.h"

mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz, bbox_t* bbox)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);

  ASSERT((bbox == NULL) || (bbox->x2 > bbox->x1));
  ASSERT((bbox == NULL) || (bbox->y2 > bbox->y1));
  ASSERT((bbox == NULL) || (bbox->z2 > bbox->z1));

  // Set up generators.
  // FIXME: Not parallel-safe!
  double dx = (bbox->x2 - bbox->x1) / nx;
  double dy = (bbox->y2 - bbox->y1) / ny;
  double dz = (bbox->z2 - bbox->z1) / nz;
  point_t* generators = malloc(sizeof(point_t) * nx * ny * nz);
  int offset = 0;
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int k = 0; k < nz; ++k, ++offset)
      {
        generators[offset].x = (i + 0.5) * dx;
        generators[offset].y = (j + 0.5) * dy;
        generators[offset].z = (k + 0.5) * dz;
      }
    }
  }

  // Create the mesh.
  mesh_t* mesh = create_voronoi_mesh_in_box(generators, nx * ny * nz,
                                            NULL, 0, bbox);

  return mesh;
}

