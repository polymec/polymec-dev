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
#include "geometry/tetrahedron.h"

static mesh_t* create_dual_mesh_from_tet_mesh(MPI_Comm comm, 
                                              mesh_t* tet_mesh,
                                              char** boundary_face_tags,
                                              int num_boundary_face_tags,
                                              char** model_edge_tags,
                                              int num_model_edge_tags,
                                              char** model_vertex_tags,
                                              int num_model_vertex_tags)
{
  // Build sets containing the indices of mesh elements identifying 
  // geometric structure (for ease of querying).
  int_unordered_set_t* boundary_tets = int_unordered_set_new();
  int_unordered_set_t* boundary_faces = int_unordered_set_new();
  for (int i = 0; i < num_boundary_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, boundary_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];
      int_unordered_set_insert(boundary_faces, face);
      int btet1 = tet_mesh->face_cells[2*face];
      int_unordered_set_insert(boundary_tets, btet1);
      int btet2 = tet_mesh->face_cells[2*face+1];
      if (btet2 != -1)
        int_unordered_set_insert(boundary_tets, btet2);
    }
  }
  int_unordered_set_t* model_edges = int_unordered_set_new();
  for (int i = 0; i < num_model_edge_tags; ++i)
  {
    int num_edges;
    int* tag = mesh_tag(tet_mesh->edge_tags, model_edge_tags[i], &num_edges);
    for (int e = 0; e < num_edges; ++e)
    {
      int edge = tag[e];
      int_unordered_set_insert(model_edges, edge);
    }
  }
  int_unordered_set_t* model_vertices = int_unordered_set_new();
  for (int i = 0; i < num_model_vertex_tags; ++i)
  {
    int num_vertices;
    int* tag = mesh_tag(tet_mesh->node_tags, model_vertex_tags[i], &num_vertices);
    for (int v = 0; v < num_vertices; ++v)
    {
      int vertex = tag[v];
      int_unordered_set_insert(model_vertices, vertex);
    }
  }

  // Allocate storage for dual vertices.
  int num_vertices = boundary_faces->size + tet_mesh->num_cells + 
                     model_edges->size + model_vertices->size;
  point_t* vertices = malloc(sizeof(point_t) * num_vertices);

  // Generate dual vertices for each of the interior tetrahedra.
  tetrahedron_t* tet = tetrahedron_new();
  for (int c = 0; c < tet_mesh->num_cells; ++c)
  {
    // The dual vertex is located at the circumcenter of the tetrahedral 
    // cell, or the point in the cell closest to it.
    point_t xc;
    tetrahedron_compute_circumcenter(tet, &xc);
    tetrahedron_compute_nearest_point(tet, &xc, &vertices[c]);
  }

  // Generate dual vertices for each of the boundary faces.
  int dv_offset = tet_mesh->num_cells;
  for (int i = 0; i < num_boundary_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, boundary_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f, ++dv_offset)
    {
      int face = tag[f];
      vertices[dv_offset] = tet_mesh->face_centers[face];
    }
  }

  // Generate a dual vertex at the midpoint of each model edge.
  for (int i = 0; i < num_model_edge_tags; ++i)
  {
    int num_edges;
    int* tag = mesh_tag(tet_mesh->edge_tags, model_edge_tags[i], &num_edges);
    for (int e = 0; e < num_edges; ++e, ++dv_offset)
    {
      int edge = tag[e];
      point_t* x1 = &tet_mesh->nodes[tet_mesh->edge_nodes[2*edge]];
      point_t* x2 = &tet_mesh->nodes[tet_mesh->edge_nodes[2*edge+1]];
      vertices[dv_offset].x = 0.5 * (x1->x + x2->x);
      vertices[dv_offset].y = 0.5 * (x1->y + x2->y);
      vertices[dv_offset].z = 0.5 * (x1->z + x2->z);
    }
  }

  // Generate a dual vertex for each model vertex.
  for (int i = 0; i < num_model_vertex_tags; ++i)
  {
    int num_vertices;
    int* tag = mesh_tag(tet_mesh->node_tags, model_vertex_tags[i], &num_vertices);
    for (int v = 0; v < num_vertices; ++v, ++dv_offset)
    {
      int vertex = tag[v];
      vertices[dv_offset] = tet_mesh->nodes[vertex];
    }
  }

  // Now that we know the various populations, build the dual mesh.
  int num_cells = 0, num_ghost_cells = 0, num_faces = 0;
  mesh_t* mesh = mesh_new(comm, num_cells, num_ghost_cells, num_faces, num_vertices);

  // Clean up.
  int_unordered_set_free(model_vertices);
  int_unordered_set_free(model_edges);
  int_unordered_set_free(boundary_tets);
  int_unordered_set_free(boundary_faces);

  return mesh;
}

mesh_t* create_dual_mesh(MPI_Comm comm, 
                         mesh_t* original_mesh,
                         char** boundary_face_tags,
                         int num_boundary_face_tags,
                         char** model_edge_tags,
                         int num_model_edge_tags,
                         char** model_vertex_tags,
                         int num_model_vertex_tags)
{
  ASSERT(num_boundary_face_tags > 0);
  ASSERT(num_boundary_face_tags > 0);

  // Currently, we only support duals of tet meshes.
  ASSERT(mesh_has_feature(original_mesh, TETRAHEDRAL));
  return create_dual_mesh_from_tet_mesh(comm, original_mesh, 
                                        boundary_face_tags, num_boundary_face_tags,
                                        model_edge_tags, num_model_edge_tags,
                                        model_vertex_tags, num_model_vertex_tags);
}

