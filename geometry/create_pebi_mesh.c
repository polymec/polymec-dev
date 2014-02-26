// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "geometry/create_pebi_mesh.h"

const char* PEBI = "perpendicular bisector";

mesh_t* create_pebi_mesh(MPI_Comm comm, 
                         point_t* cell_centers, real_t* cell_volumes, int num_cells,
                         int* faces, real_t* face_areas, point_t* face_centers, 
                         int num_faces)
{
  // Check input.
  ASSERT(cell_centers != NULL);
  ASSERT(cell_volumes != NULL);
  ASSERT(num_cells > 0);
  ASSERT(faces != NULL);
  ASSERT(face_areas != NULL);
  ASSERT(num_faces >= 0);
#ifndef NDEBUG
  for (int f = 0; f < num_faces; ++f)
  {
    ASSERT(faces[2*f] >= 0);
    ASSERT((faces[2*f+1] >= 0) || (faces[2*f+1] == -1));
  }
#endif

  mesh_t* mesh = mesh_new(comm, num_cells, 0, num_faces, 0, 0);

  // Transcribe the mesh cell centers, which are the only connection to 
  // spatial geometry.
  memcpy(mesh->cell_centers, cell_centers, sizeof(point_t)*num_cells);

  // Copy over the face-cell connectivity directly.
  memcpy(mesh->face_cells, faces, 2*sizeof(int)*num_faces);

  // Copy over the face areas directly.
  memcpy(mesh->face_areas, face_areas, sizeof(real_t)*num_faces);

  // Go through the list of faces and count the faces attached to each cell,
  // storing the tally in the set of cell face offsets.
  for (int f = 0; f < num_faces; ++f)
  {
    mesh->cell_face_offsets[faces[2*f]] += 1;
    if (faces[2*f+1] != -1)
      mesh->cell_face_offsets[faces[2*f]] += 1;
  }
  // Convert these face offsets to compressed row storage format.
  for (int c = 1; c <= num_cells; ++c)
    mesh->cell_face_offsets[c] += mesh->cell_face_offsets[c-1];

  // Now fill the mesh's cell_faces array.
  int* cell_face_count = malloc(sizeof(int) * num_cells);
  memset(cell_face_count, 0, sizeof(int) * num_cells);
  mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int) * mesh->cell_face_offsets[num_cells], 0);
  for (int f = 0; f < num_faces; ++f)
  {
    int c1 = faces[2*f];
    mesh->cell_faces[c1] = f;
    ++cell_face_count[c1];
    int c2 = faces[2*f+1];
    if (c2 != -1)
    {
      mesh->cell_faces[c2] = f;
      ++cell_face_count[c2];
    }
  }

  // Set the cell volumes.
  memcpy(mesh->cell_volumes, cell_volumes, sizeof(real_t)*num_cells);

  // Set or compute face centers.
  // We compute information for interior faces first.
  for (int f = 0; f < num_faces; ++f)
  {
    point_t* xf = &mesh->face_centers[f];
    vector_t* nf = &mesh->face_normals[f];
    int c1 = mesh->face_cells[2*f];
    int c2 = mesh->face_cells[2*f];
    if (c2 != -1) // Interior face
    {
      point_t* xc1 = &mesh->cell_centers[c1];
      point_t* xc2 = &mesh->cell_centers[c2];
      if (face_centers == NULL)
      {
        // Assume each face center lies at the midpoint between its cells.
        xf->x = 0.5 * (xc1->x + xc2->x);
        xf->y = 0.5 * (xc1->y + xc2->y);
        xf->z = 0.5 * (xc1->z + xc2->z);
      }
      else
      {
        xf->x = face_centers[f].x;
        xf->y = face_centers[f].y;
        xf->z = face_centers[f].z;
      }

      // The face normal should connect xc1 and xc2.
      point_displacement(xc1, xc2, nf);
      vector_normalize(nf);
    }

    // Now use the existing information to compute information for 
    // boundary faces.
    for (int f = 0; f < num_faces; ++f)
    {
      int c1 = mesh->face_cells[2*f];
      int c2 = mesh->face_cells[2*f];
      if (c2 == -1)
      {
        // Estimate the cell-face distance by assuming an isotropic cell.
        real_t V = mesh->cell_volumes[c1];
        real_t d = pow(V, 1.0/3.0);

        // Form the normal vector for the face by assuming that all face 
        // normals sum to zero.
        vector_t* nf = &mesh->face_normals[f];
        nf->x = nf->y = nf->z = 0.0;
        for (int ff = mesh->cell_face_offsets[c1]; ff < mesh->cell_face_offsets[c1+1]; ++ff)
        {
          vector_t* nff = &mesh->face_normals[ff];
          nf->x -= nff->x;
          nf->y -= nff->y;
          nf->z -= nff->z;
        }

        // Compute the face center.
        if (face_centers == NULL)
        {
          point_t* xc = &mesh->cell_centers[c1];
          point_t* xf = &mesh->face_centers[f];
          xf->x = xc->x + d*nf->x;
          xf->y = xc->y + d*nf->y;
          xf->z = xc->z + d*nf->z;
        }
        else 
        {
          xf->x = face_centers[f].x;
          xf->y = face_centers[f].y;
          xf->z = face_centers[f].z;
        }
      }
    }
  }

  mesh_add_feature(mesh, PEBI);
  return mesh;
}

mesh_t* create_pebi_mesh_from_unstructured_mesh(mesh_t* mesh)
{
  return create_pebi_mesh(mesh->comm, mesh->cell_centers, mesh->cell_volumes, 
                          mesh->num_cells, mesh->face_cells, mesh->face_areas, 
                          mesh->face_centers, mesh->num_faces);
}
