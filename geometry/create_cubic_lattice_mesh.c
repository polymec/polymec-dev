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
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/create_rectilinear_lattice_mesh.h"
#include "geometry/cubic_lattice.h"

mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz, bbox_t* bbox)
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
  double x1 = bbox->x1;
  double y1 = bbox->y1;
  double z1 = bbox->z1;
  double dx = Lx/nx, dy = Ly/ny, dz = Lz/nz;

  // Create a uniform rectilinear lattice mesh!
  double* xs = malloc(sizeof(double) * (nx+1));
  double* ys = malloc(sizeof(double) * (ny+1));
  double* zs = malloc(sizeof(double) * (nz+1));
  for (int i = 0; i <= nx; ++i)
    xs[i] = i*dx;
  for (int i = 0; i <= ny; ++i)
    ys[i] = i*dy;
  for (int i = 0; i <= nz; ++i)
    zs[i] = i*dz;

  mesh_t* mesh = create_rectilinear_lattice_mesh(xs, nx+1, ys, ny+1, zs, nz+1);

  // Clean up.
  free(zs);
  free(ys);
  free(xs);

  return mesh;
}

void tag_cubic_lattice_mesh_faces(mesh_t* mesh, 
                                  int nx, int ny, int nz,
                                  const char* x1_tag, 
                                  const char* x2_tag, 
                                  const char* y1_tag,
                                  const char* y2_tag,
                                  const char* z1_tag,
                                  const char* z2_tag)
{
  // Tag the boundaries of the mesh.
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int* x1tag = mesh_create_tag(mesh->face_tags, x1_tag, ny*nz);
  int* x2tag = mesh_create_tag(mesh->face_tags, x2_tag, ny*nz);
  for (int j = 0; j < ny; ++j)
  {
    for (int k = 0; k < nz; ++k)
    {
      x1tag[nz*j + k] = cubic_lattice_x_face(lattice, 0, j, k);
      x2tag[nz*j + k] = cubic_lattice_x_face(lattice, nx, j, k);
    }
  }

  int* y1tag = mesh_create_tag(mesh->face_tags, y1_tag, nx*nz);
  int* y2tag = mesh_create_tag(mesh->face_tags, y2_tag, nx*nz);
  for (int i = 0; i < nx; ++i)
  {
    for (int k = 0; k < nz; ++k)
    {
      y1tag[nz*i + k] = cubic_lattice_y_face(lattice, i, 0, k);
      y2tag[nz*i + k] = cubic_lattice_y_face(lattice, i, ny, k);
    }
  }

  int* z1tag = mesh_create_tag(mesh->face_tags, z1_tag, nx*ny);
  int* z2tag = mesh_create_tag(mesh->face_tags, z2_tag, nx*ny);
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      z1tag[ny*i + j] = cubic_lattice_z_face(lattice, i, j, 0);
      z2tag[ny*i + j] = cubic_lattice_z_face(lattice, i, j, nz);
    }
  }
  lattice = NULL;
}

