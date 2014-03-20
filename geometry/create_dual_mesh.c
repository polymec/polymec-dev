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
#include "core/unordered_map.h"
#include "geometry/tetrahedron.h"
#include "geometry/plane.h"
#include "geometry/polygon.h"
#include "geometry/create_dual_mesh.h"

// This type and its comparator (below) are used with qsort to order face nodes.
typedef struct 
{
  double angle;
  int index;
} dual_face_node_ordering_t;

static inline int dual_face_node_cmp(const void* l, const void* r)
{
  const dual_face_node_ordering_t* fl = l;
  const dual_face_node_ordering_t* fr = r;
  return (fl->angle < fr->angle) ? -1 :
         (fl->angle > fr->angle) ?  1 : 0;             
}


// Given a set of nodes in a plane (with two identified as "endpoints", 
// order them from the first endpoint to the second, placing their indices 
// (corresponding to their locations in dual_nodes) in dual_node_indices.
static void order_nodes_of_dual_face(sp_func_t* plane, 
                                     int* endpoint_indices, 
                                     point_t* dual_nodes, 
                                     int num_nodes, 
                                     int* dual_node_indices)
{
  ASSERT(num_nodes >= 2);
  ASSERT(endpoint_indices[0] != -1);
  ASSERT(endpoint_indices[1] != -1);

  if (num_nodes == 2)
  {
    // If there are only two nodes, we're done.
    dual_node_indices[0] = endpoint_indices[0];
    dual_node_indices[1] = endpoint_indices[1];
  }
  else if (num_nodes == 3)
  {
    // If there are only three nodes, there's only one in between them.
    dual_node_indices[0] = endpoint_indices[0];
    dual_node_indices[2] = endpoint_indices[1];
    for (int i = 0; i < 3; ++i)
    {
      if ((i != dual_node_indices[0]) && (i != dual_node_indices[2]))
      {
        dual_node_indices[1] = i;
        break;
      }
    }
  }

  // Otherwise, we have to sort the remaining nodes.
  dual_face_node_ordering_t tweener_nodes[num_nodes-2];
  int j = 0;
  for (int i = 0; i < num_nodes; ++i)
  {
    if ((i != dual_node_indices[0]) && (i != dual_node_indices[2]))
    {
      // Project these nodes to our plane.
      point2_t xi;
      plane_project(plane, &dual_nodes[i], &xi);
      tweener_nodes[j].angle = atan2(xi.y, xi.x);
      tweener_nodes[j].index = i;
    }
  }
  qsort(tweener_nodes, (size_t)(num_nodes-2), sizeof(dual_face_node_ordering_t), dual_face_node_cmp);

  // Put everything back into place.
  dual_node_indices[0] = endpoint_indices[0];
  for (int i = 0; i < num_nodes-2; ++i)
    dual_node_indices[i+1] = tweener_nodes[i].index;
  dual_node_indices[num_nodes-1] = endpoint_indices[0];
}

static mesh_t* create_dual_mesh_from_tet_mesh(MPI_Comm comm, 
                                              mesh_t* tet_mesh,
                                              char** external_model_face_tags,
                                              int num_external_model_face_tags,
                                              char** internal_model_face_tags,
                                              int num_internal_model_face_tags,
                                              char** model_edge_tags,
                                              int num_model_edge_tags,
                                              char** model_vertex_tags,
                                              int num_model_vertex_tags)
{
  // Build sets containing the indices of mesh elements identifying 
  // geometric structure (for ease of querying).
  int_unordered_set_t* external_boundary_tets = int_unordered_set_new();
  int_unordered_set_t* external_model_faces = int_unordered_set_new();
  int_unordered_set_t* external_model_face_edges = int_unordered_set_new();
  for (int i = 0; i < num_external_model_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, external_model_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];
      int_unordered_set_insert(external_model_faces, face);

      int btet1 = tet_mesh->face_cells[2*face];
      int_unordered_set_insert(external_boundary_tets, btet1);
      int btet2 = tet_mesh->face_cells[2*face+1];
      if (btet2 != -1)
        int_unordered_set_insert(external_boundary_tets, btet2);

      int pos = 0, edge;
      while (mesh_face_next_edge(tet_mesh, face, &pos, &edge))
        int_unordered_set_insert(external_model_face_edges, edge);
    }
  }
  int_unordered_set_t* internal_boundary_tets = int_unordered_set_new();
  int_unordered_set_t* internal_model_faces = int_unordered_set_new();
  int_unordered_set_t* internal_model_face_edges = int_unordered_set_new();
  for (int i = 0; i < num_internal_model_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, internal_model_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];
      int_unordered_set_insert(internal_model_faces, face);

      int btet1 = tet_mesh->face_cells[2*face];
      int_unordered_set_insert(internal_boundary_tets, btet1);
      int btet2 = tet_mesh->face_cells[2*face+1];
      if (btet2 != -1)
        int_unordered_set_insert(internal_boundary_tets, btet2);

      int pos = 0, edge;
      while (mesh_face_next_edge(tet_mesh, face, &pos, &edge))
        int_unordered_set_insert(internal_model_face_edges, edge);
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
  int num_dual_nodes = external_model_faces->size + 
                       internal_model_faces->size + tet_mesh->num_cells + 
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

  // Generate dual vertices for each of the model faces.
  for (int i = 0; i < num_external_model_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, external_model_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f, ++dv_offset)
    {
      int face = tag[f];
      dual_mesh->nodes[dv_offset] = tet_mesh->face_centers[face];
    }
  }
  for (int i = 0; i < num_internal_model_face_tags; ++i)
  {
    int num_faces;
    int* tag = mesh_tag(tet_mesh->face_tags, internal_model_face_tags[i], &num_faces);
    for (int f = 0; f < num_faces; ++f, ++dv_offset)
    {
      int face = tag[f];
      dual_mesh->nodes[dv_offset] = tet_mesh->face_centers[face];
    }
  }

  // Generate a dual vertex at the midpoint of each model edge.
  int_int_unordered_map_t* dual_node_for_edge = int_int_unordered_map_new();
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
      int_int_unordered_map_insert(dual_node_for_edge, e, dv_offset);
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
    int_unordered_set_t* cells_for_edge = primal_cells_for_edge[edge];
    ASSERT(cells_for_edge != NULL);

    // Is this edge a model edge?
    bool is_external_face_edge = int_unordered_set_contains(external_model_face_edges, edge);
    bool is_internal_face_edge = int_unordered_set_contains(internal_model_face_edges, edge);
    bool is_model_edge = int_unordered_set_contains(model_edges, edge);

    if (is_external_face_edge || is_internal_face_edge)
    {
      if (is_external_face_edge)
      {
        // This primal edge belongs to an external model face,
        // so it lies on the outside of the domain. The corresponding dual 
        // face is bounded by dual nodes created from the primal cells 
        // bounding the edge. We want to order these dual nodes starting at 
        // one boundary cell and finishing at the other. So we extract the 
        // indices of the dual nodes and then pick out the endpoints.
        int num_nodes = cells_for_edge->size;
        int pos = 0, cell, c = 0;
        point_t dual_nodes[num_nodes];
        int dual_node_indices[num_nodes];
        int endpoint_indices[] = {-1, -1};
        while (int_unordered_set_next(cells_for_edge, &pos, &cell))
        {
          dual_nodes[c] = tet_mesh->cell_centers[cell];
          dual_node_indices[c] = cell;
          if (int_unordered_set_contains(external_boundary_tets, cell))
          {
            if (endpoint_indices[0] == -1)
              endpoint_indices[0] = c;
            else
              endpoint_indices[1] = c;
          }
          ++c;
        }

        // Find a vector connecting the nodes of this edge. This orients 
        // the face.
        vector_t edge_vector;
        point_t* x1 = &tet_mesh->nodes[tet_mesh->edge_nodes[2*edge]];
        point_t* x2 = &tet_mesh->nodes[tet_mesh->edge_nodes[2*edge+1]];
        point_displacement(x1, x2, &edge_vector);

        // Order the nodes of this dual face.
        sp_func_t* edge_plane = plane_new(&edge_vector, x1);
        order_nodes_of_dual_face(edge_plane, endpoint_indices, dual_nodes, 
            num_nodes, dual_node_indices);

        // Update the dual mesh's face->node connectivity metadata.
        int num_face_nodes = (is_model_edge) ? num_nodes + 1 : num_nodes;
        dual_mesh->face_node_offsets[df_offset+1] = dual_mesh->face_node_offsets[df_offset] + num_face_nodes;
        int_array_t* face_nodes = int_array_new();
        nodes_for_dual_face[df_offset] = face_nodes;
        int_array_resize(face_nodes, num_face_nodes);
        memcpy(face_nodes->data, dual_node_indices, sizeof(int)*num_nodes);

        // If the edge is a model edge, stick the primal edge's node at the end 
        // of the list of dual face nodes.
        if (is_model_edge)
        {
          int* dual_node_from_edge_p = int_int_unordered_map_get(dual_node_for_edge, edge);
          ASSERT(dual_node_from_edge_p != NULL);
          int dual_node_from_edge = *dual_node_from_edge_p;
          face_nodes->data[num_nodes] = dual_node_from_edge;
        }
        ++df_offset;
      }
      else
      {
        // This primal edge belongs to an internal model face,
        // so it lies on an interface between two regions within the domain. 
        // We create two dual faces for this edge (one for each region), 
        // using a procedure very similar to the one we used for external 
        // edges above.

        // Dump the IDs of the cells attached to the edge into a single array.
        int num_cells = cells_for_edge->size;
        int pos = 0, cell, c = 0;
        point_t dual_nodes[num_cells];
        int dual_node_indices[num_cells];
        while (int_unordered_set_next(cells_for_edge, &pos, &cell))
        {
          dual_nodes[c] = tet_mesh->cell_centers[cell];
          dual_node_indices[c] = cell;
          ++c;
        }

        // Since this is an internal interface edge, the dual nodes 
        // corresponding to these cells form a polygon around the edge. We can 
        // arrange the nodes for the two faces (stuck together) into a polygon 
        // using the "star" algorithm and then retrieve them (in order) from 
        // the polygon.
        polygon_t* dual_polygon = polygon_giftwrap(dual_nodes, num_cells);
//        polygon_t* dual_polygon = polygon_star(dual_nodes, num_cells);
        int* ordering = polygon_ordering(dual_polygon);

        // Now we just need to apportion the right nodes to the right faces.
        int start_index1 = -1, start_index2 = -1, stop_index1 = -1, stop_index2 = -1;
        for (int i = 0; i < num_cells; ++i)
        {
          int this_cell = dual_node_indices[ordering[i]];
          int next_cell = dual_node_indices[ordering[(i+1)%num_cells]];
          // If this_cell and next_cell share a face that is an internal 
          // model face, they are on the opposite side of the interface.
          if (int_unordered_set_contains(internal_boundary_tets, this_cell))
          {
          }
        }

        // Update the dual mesh's face->node connectivity metadata for 
        // both faces.
        int num_nodes1 = stop_index1 - start_index1 + 1;
        int num_face1_nodes = (is_model_edge) ? num_nodes1 + 1 : num_nodes1;
        dual_mesh->face_node_offsets[df_offset+1] = dual_mesh->face_node_offsets[df_offset] + num_face1_nodes;
        int_array_t* face1_nodes = int_array_new();
        nodes_for_dual_face[df_offset] = face1_nodes;
        int_array_resize(face1_nodes, num_nodes1);
        for (int i = start_index1; i <= stop_index1; ++i)
        {
          int j = (start_index1 + i) % num_cells;
          face1_nodes->data[i] = dual_node_indices[ordering[j]];
        }

        int num_nodes2 = stop_index2 - start_index2 + 1;
        int num_face2_nodes = (is_model_edge) ? num_nodes2 + 1 : num_nodes2;
        dual_mesh->face_node_offsets[df_offset+2] = dual_mesh->face_node_offsets[df_offset+1] + num_face2_nodes;
        int_array_t* face2_nodes = int_array_new();
        nodes_for_dual_face[df_offset+1] = face2_nodes;
        int_array_resize(face2_nodes, num_nodes2);
        for (int i = 0; i <= num_nodes2; ++i)
        {
          int j = (start_index2 + i) % num_cells;
          face1_nodes->data[i] = dual_node_indices[ordering[j]];
        }

        // If the edge is a model edge, stick the primal edge's node at the end 
        // of each of the lists of dual face nodes.
        if (is_model_edge)
        {
          int* dual_node_from_edge_p = int_int_unordered_map_get(dual_node_for_edge, edge);
          ASSERT(dual_node_from_edge_p != NULL);
          int dual_node_from_edge = *dual_node_from_edge_p;
          face1_nodes->data[num_nodes1] = dual_node_from_edge;
          face2_nodes->data[num_nodes2] = dual_node_from_edge;
        }
        df_offset += 2;
      }
    }
    else
    {
      // This edge is on the interior of the domain, so it is only bounded 
      // by cells.

      // Dump the cell centers into an array.
      int num_cells = cells_for_edge->size;
      int pos = 0, cell, c = 0;
      point_t dual_nodes[num_cells];
      while (int_unordered_set_next(cells_for_edge, &pos, &cell))
        dual_nodes[c++] = tet_mesh->cell_centers[cell];

      // Update the dual mesh's connectivity metadata.
      dual_mesh->face_node_offsets[df_offset+1] = dual_mesh->face_node_offsets[df_offset] + num_cells;

      // Since this is an interior edge, the dual nodes corresponding to 
      // these cells form a convex polygon around the edge. We can arrange 
      // the nodes into a convex polygon using the gift-wrapping algorithm.
      polygon_t* dual_polygon = polygon_giftwrap(dual_nodes, num_cells);
      int_array_t* face_nodes = int_array_new();
      int_array_resize(face_nodes, num_cells);
      memcpy(face_nodes->data, polygon_ordering(dual_polygon), sizeof(int)*num_cells);
      nodes_for_dual_face[df_offset] = face_nodes;
      dual_polygon = NULL;
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
  int_int_unordered_map_free(dual_node_for_edge);
  int_unordered_set_free(model_vertices);
  int_unordered_set_free(model_edges);
  int_unordered_set_free(external_boundary_tets);
  int_unordered_set_free(external_model_faces);
  int_unordered_set_free(external_model_face_edges);
  int_unordered_set_free(internal_boundary_tets);
  int_unordered_set_free(internal_model_faces);
  int_unordered_set_free(internal_model_face_edges);

  // Compute mesh geometry.
  mesh_compute_geometry(dual_mesh);

  return dual_mesh;
}

mesh_t* create_dual_mesh(MPI_Comm comm, 
                         mesh_t* original_mesh,
                         char** external_model_face_tags,
                         int num_external_model_face_tags,
                         char** internal_model_face_tags,
                         int num_internal_model_face_tags,
                         char** model_edge_tags,
                         int num_model_edge_tags,
                         char** model_vertex_tags,
                         int num_model_vertex_tags)
{
  ASSERT(num_external_model_face_tags > 0);
  ASSERT(num_model_edge_tags > 0);
  ASSERT(num_model_vertex_tags > 0);

  // Currently, we only support duals of tet meshes.
  ASSERT(mesh_has_feature(original_mesh, TETRAHEDRAL));
  return create_dual_mesh_from_tet_mesh(comm, original_mesh, 
                                        external_model_face_tags, num_external_model_face_tags,
                                        internal_model_face_tags, num_internal_model_face_tags,
                                        model_edge_tags, num_model_edge_tags,
                                        model_vertex_tags, num_model_vertex_tags);
}

