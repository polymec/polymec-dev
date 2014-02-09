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
                         point_t* cell_centers, int num_cells,
                         int* faces, real_t* face_areas, int num_faces)
{
#ifndef NDEBUG
  // Check input.
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

  return mesh;
}

