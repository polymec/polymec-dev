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
#include "geometry/tetrahedron.h"
#include "geometry/polygon.h"
#include "geometry/create_dual_mesh.h"

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

  // Count up the dual mesh entities.
  int num_dual_nodes = boundary_faces->size + tet_mesh->num_cells + 
                       model_edges->size + model_vertices->size;
  int num_dual_faces = tet_mesh->num_edges - model_edges->size;
  int num_dual_cells = 0, num_dual_ghost_cells = 0; 
  // FIXME

  // Now that we know the various populations, build the dual mesh.
  mesh_t* dual_mesh = mesh_new(comm, num_dual_cells, num_dual_ghost_cells, 
                               num_dual_faces, num_dual_nodes);

  // Generate dual vertices for each of the interior tetrahedra.
  tetrahedron_t* tet = tetrahedron_new();
  int dv_offset = 0;
  for (int c = 0; c < tet_mesh->num_cells; ++c, ++dv_offset)
  {
    // The dual vertex is located at the circumcenter of the tetrahedral 
    // cell, or the point in the cell closest to it.
    point_t xc;
    tetrahedron_compute_circumcenter(tet, &xc);
    tetrahedron_compute_nearest_point(tet, &xc, &dual_mesh->nodes[dv_offset]);
  }

  // Generate dual vertices for each of the boundary faces.
  for (int i = 0; i < num_boundary_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, boundary_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f, ++dv_offset)
    {
      int face = tag[f];
      dual_mesh->nodes[dv_offset] = tet_mesh->face_centers[face];
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
      point_t* n = &dual_mesh->nodes[dv_offset];
      n->x = 0.5 * (x1->x + x2->x);
      n->y = 0.5 * (x1->y + x2->y);
      n->z = 0.5 * (x1->z + x2->z);
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
      dual_mesh->nodes[dv_offset] = tet_mesh->nodes[vertex];
    }
  }

  // Now generate dual faces corresponding to primal edges. Each primal edge 
  // is surrounded by primal cells that correspond to dual vertices, so we 
  // have to build the edge->cell connectivity and then make sure that the 
  // cells around an edge are ordered in a well-defined manner.
  int_unordered_set_t** primal_cells_for_edge = malloc(sizeof(int_unordered_set_t*) * tet_mesh->num_edges);
  memset(primal_cells_for_edge, 0, sizeof(int_unordered_set_t*) * tet_mesh->num_edges);
  for (int cell = 0; cell < tet_mesh->num_cells; ++cell)
  {
    // Loop over faces in the cell and associate this cell with each of the edges.
    int pos = 0, face;
    while (mesh_cell_next_face(tet_mesh, cell, &pos, &face))
    {
      int pos1 = 0, edge;
      while (mesh_face_next_edge(tet_mesh, face, &pos1, &edge))
      {
        int_unordered_set_t* cells_for_edge = primal_cells_for_edge[edge];
        if (cells_for_edge == NULL)
        {
          cells_for_edge = int_unordered_set_new();
          primal_cells_for_edge[edge] = cells_for_edge;
        }
        int_unordered_set_insert(cells_for_edge, cell);
      }
    }
  }

  // Now step over the primal edges and create dual faces.
  int df_offset = 0;
  dual_mesh->face_node_offsets[0] = 0;
  int_array_t** nodes_for_dual_face = malloc(sizeof(int_array_t*) * num_dual_faces);
  memset(nodes_for_dual_face, 0, sizeof(int_array_t*) * num_dual_faces);
  for (int edge = 0; edge < tet_mesh->num_edges; ++edge)
  {
    if (int_unordered_set_contains(model_edges, edge))
    {
      // This edge is a model edge, so it lies on the exterior of the 
      // domain. We have to be careful.
    }
    else
    {
      int_unordered_set_t* cells_for_edge = primal_cells_for_edge[edge];
      ASSERT(cells_for_edge != NULL);

      // Dump the cell IDs into an array.
      int num_cells = cells_for_edge->size;
      int pos = 0, cell, c = 0, primal_cells[num_cells];
      point_t dual_vertices[num_cells];
      while (int_unordered_set_next(cells_for_edge, &pos, &cell))
      {
        dual_vertices[c] = tet_mesh->cell_centers[cell];
        primal_cells[c] = cell;
        ++c;
      }

      // Update the dual mesh's connectivity metadata.
      dual_mesh->face_node_offsets[df_offset+1] = dual_mesh->face_node_offsets[df_offset] + num_cells;

      // Since this is an interior edge, the dual vertices corresponding to 
      // these cells form a convex polygon around the edge. We can arrange 
      // the vertices into a convex polygon using the gift-wrapping algorithm.
      polygon_t* dual_polygon = polygon_giftwrap(dual_vertices, num_cells);
      int_array_t* face_nodes = int_array_new();
      int_array_resize(face_nodes, num_cells);
      memcpy(face_nodes->data, polygon_ordering(dual_polygon), sizeof(int)*num_cells);
      nodes_for_dual_face[df_offset] = face_nodes;
      ++df_offset;
    }
  }

  // Create dual cells.
  int dc_offset = 0;
  int_array_t** faces_for_dual_cell = malloc(sizeof(int_array_t*) * num_dual_cells);
  memset(faces_for_dual_cell, 0, sizeof(int_array_t*) * num_dual_cells);
  // FIXME

  // Allocate mesh connectivity storage and move all the data into place.
  mesh_reserve_connectivity_storage(dual_mesh);
  for (int c = 0; c < num_dual_cells; ++c)
  {
    int_array_t* cell_faces = faces_for_dual_cell[c];
    memcpy(&dual_mesh->cell_faces[dual_mesh->cell_face_offsets[c]], cell_faces->data, sizeof(int)*cell_faces->size);
    for (int f = 0; f < cell_faces->size; ++f)
    {
      int face = cell_faces->data[f];
      if (dual_mesh->face_cells[2*face] == -1)
        dual_mesh->face_cells[2*face] = c;
      else
        dual_mesh->face_cells[2*face+1] = c;
    }
  }
  for (int f = 0; f < num_dual_faces; ++f)
  {
    int_array_t* face_nodes = nodes_for_dual_face[f];
    memcpy(&dual_mesh->face_nodes[dual_mesh->face_node_offsets[f]], face_nodes->data, sizeof(int)*face_nodes->size);
  }

  // Clean up.
  for (int c = 0; c < num_dual_cells; ++c)
    int_array_free(faces_for_dual_cell[c]);
  free(faces_for_dual_cell);
  for (int f = 0; f < num_dual_faces; ++f)
    int_array_free(nodes_for_dual_face[f]);
  free(nodes_for_dual_face);
  for (int e = 0; e < tet_mesh->num_edges; ++e)
    int_unordered_set_free(primal_cells_for_edge[e]);
  free(primal_cells_for_edge);
  int_unordered_set_free(model_vertices);
  int_unordered_set_free(model_edges);
  int_unordered_set_free(boundary_tets);
  int_unordered_set_free(boundary_faces);

  // Compute mesh geometry.
  mesh_compute_geometry(dual_mesh);

  return dual_mesh;
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

