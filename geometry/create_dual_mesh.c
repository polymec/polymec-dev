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

#include "core/unordered_set.h"
#include "geometry/create_dual_mesh.h"

static mesh_t* create_dual_mesh_from_tet_mesh(MPI_Comm comm, 
                                              mesh_t* tet_mesh,
                                              char** boundary_face_tags,
                                              int num_boundary_face_tags)
{
  // Build sets containing the indices of boundary tetrahedra and faces 
  // (for ease of querying).
  int_unordered_set_t* boundary_tets = int_unordered_set_new();
  int_unordered_set_t* boundary_faces = int_unordered_set_new();
  for (int i = 0; i < num_boundary_face_tags; ++i)
  {
    int num_bfaces;
    int* btag = mesh_tag(tet_mesh->face_tags, boundary_face_tags[i], &num_bfaces);
    for (int f = 0; f < num_bfaces; ++f)
    {
      int bface = btag[f];
      int_unordered_set_insert(boundary_faces, bface);
      int btet1 = tet_mesh->face_cells[2*bface];
      int_unordered_set_insert(boundary_tets, btet1);
      int btet2 = tet_mesh->face_cells[2*bface+1];
      if (btet2 != -1)
        int_unordered_set_insert(boundary_tets, btet2);
    }
  }

  // Allocate storage for dual vertices.
  int num_vertices = boundary_faces->size + tet_mesh->num_cells;
  point_t* vertices = malloc(sizeof(point_t) * num_vertices);

  // Generate dual vertices for each of the interior tetrahedra.
  for (int c = 0; c < tet_mesh->num_cells; ++c)
  {
    if (int_unordered_set_contains(boundary_tets, c))
    {
      // This tet is on the boundary and is probably not well-centered.
    }
    else
    {
      // This is an interior tet and is well-centered.
    }
  }

  // Generate dual vertices for each of the boundary faces.
  for (int i = 0; i < num_boundary_face_tags; ++i)
  {
    int num_bfaces;
    int* btag = mesh_tag(tet_mesh->face_tags, boundary_face_tags[i], &num_bfaces);
    for (int f = 0; f < num_bfaces; ++f)
    {
      int bface = btag[f];
    }
  }

  // Now that we know the various populations, build the dual mesh.
  int num_cells = 0, num_ghost_cells = 0, num_faces = 0, num_edges = 0;
  mesh_t* mesh = mesh_new(comm, num_cells, num_ghost_cells, num_faces,
                          num_edges, num_vertices);

  // Clean up.
  int_unordered_set_free(boundary_tets);
  int_unordered_set_free(boundary_faces);

  return mesh;
}

mesh_t* create_dual_mesh(MPI_Comm comm, 
                         mesh_t* original_mesh,
                         char** boundary_face_tags,
                         int num_boundary_face_tags)
{
  ASSERT(num_boundary_face_tags > 0);
  ASSERT(num_boundary_face_tags > 0);

  // Currently, we only support duals of tet meshes.
  ASSERT(mesh_has_feature(original_mesh, TETRAHEDRAL));
  return create_dual_mesh_from_tet_mesh(comm, original_mesh, boundary_face_tags, num_boundary_face_tags);
}

