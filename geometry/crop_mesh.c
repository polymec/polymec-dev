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
#include "geometry/crop_mesh.h"

mesh_t* crop_mesh(mesh_t* mesh, sp_func_t* boundary_func)
{
  // Mark the cells outside the boundary.
  int_unordered_set_t* outside_cells = int_unordered_set_new();
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    real_t dist;
    sp_func_eval(boundary_func, &mesh->cell_centers[cell], &dist);
    if (dist > 0.0)
      int_unordered_set_insert(outside_cells, cell);
  }

  // Count up the remaining faces and nodes, and make a list of boundary faces.
  int_unordered_set_t* remaining_ghosts = int_unordered_set_new();
  int_unordered_set_t* boundary_faces = int_unordered_set_new();
  int_int_unordered_map_t* remaining_nodes = int_int_unordered_map_new();
  int_int_unordered_map_t* remaining_faces = int_int_unordered_map_new();
  int new_face_index = 0, new_node_index = 0;
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    if (int_unordered_set_contains(outside_cells, cell)) continue;

    int pos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      if (!int_int_unordered_map_contains(remaining_faces, face))
        int_int_unordered_map_insert(remaining_faces, face, new_face_index++);
      int opp_cell = mesh_face_opp_cell(mesh, face, cell);
      if ((opp_cell == -1) || int_unordered_set_contains(outside_cells, opp_cell))
        int_unordered_set_insert(boundary_faces, face);
      else if (opp_cell >= mesh->num_cells)
        int_unordered_set_insert(remaining_ghosts, opp_cell);
      int pos1 = 0, node;
      while (mesh_face_next_node(mesh, face, &pos1, &node))
      {
        if (!int_int_unordered_map_contains(remaining_nodes, node))
          int_int_unordered_map_insert(remaining_nodes, node, new_node_index++);
      }
    }
  }

  // Convert the face and node maps to arrays.
  int* face_map = polymec_malloc(sizeof(int) * mesh->num_faces);
  for (int i = 0; i < mesh->num_faces; ++i)
    face_map[i] = -1;
  {
    int pos = 0, old_face, new_face;
    while (int_int_unordered_map_next(remaining_faces, &pos, &old_face, &new_face))
      face_map[old_face] = new_face;
  }
  int* node_map = polymec_malloc(sizeof(int) * mesh->num_nodes);
  for (int i = 0; i < mesh->num_nodes; ++i)
    node_map[i] = -1;
  {
    int pos = 0, old_node, new_node;
    while (int_int_unordered_map_next(remaining_nodes, &pos, &old_node, &new_node))
      node_map[old_node] = new_node;
  }

  // Do a partial Clean up.
  int num_new_faces = remaining_faces->size;
  int num_new_nodes = remaining_nodes->size;
  int num_new_ghosts = remaining_ghosts->size;
  int_int_unordered_map_free(remaining_nodes);
  int_int_unordered_map_free(remaining_faces);
  int_unordered_set_free(remaining_ghosts);

  // Make a new mesh.
  mesh_t* cropped_mesh = mesh_new(mesh->comm, 
                                  mesh->num_cells - outside_cells->size,
                                  num_new_ghosts, num_new_faces, num_new_nodes);

  // Reserve storage.
  int new_cell = 0;
  cropped_mesh->cell_face_offsets[0] = 0;
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    if (int_unordered_set_contains(outside_cells, cell)) continue;
    int num_cell_faces = mesh_cell_num_faces(mesh, cell);
    cropped_mesh->cell_face_offsets[new_cell+1] = cropped_mesh->cell_face_offsets[new_cell] + num_cell_faces;
    ++new_cell;
  }
  cropped_mesh->face_node_offsets[0] = 0;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int new_face = face_map[f];
    if (new_face != -1)
    {
      int num_face_nodes = mesh_face_num_nodes(mesh, f);
      cropped_mesh->face_node_offsets[new_face+1] = num_face_nodes;
    }
  }
  for (int f = 0; f < cropped_mesh->num_faces; ++f)
    cropped_mesh->face_node_offsets[f+1] += cropped_mesh->face_node_offsets[f];
  mesh_reserve_connectivity_storage(cropped_mesh);

  // Hook up the faces and cells.
  new_cell = 0;
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    if (int_unordered_set_contains(outside_cells, cell)) continue;
    int num_cell_faces = mesh_cell_num_faces(mesh, cell);
    for (int f = 0; f < num_cell_faces; ++f)
    {
      int old_face = mesh->cell_faces[mesh->cell_face_offsets[cell]+f];
      int pos_old_face = (old_face >= 0) ? old_face : ~old_face;
      int new_face = face_map[pos_old_face];
      if (cropped_mesh->face_cells[2*new_face] == -1)
        cropped_mesh->face_cells[2*new_face] = new_cell;
      else
        cropped_mesh->face_cells[2*new_face+1] = new_cell;
      if (old_face < 0) new_face = ~new_face;
      cropped_mesh->cell_faces[cropped_mesh->cell_face_offsets[new_cell]+f] = new_face;
    }
    ++new_cell;
  }
  ASSERT(new_cell == cropped_mesh->num_cells);

  // Hook up faces and nodes.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int new_face = face_map[f];
    if (new_face != -1)
    {
      int num_face_nodes = mesh_face_num_nodes(mesh, f);
      for (int n = 0; n < num_face_nodes; ++n)
      {
        int old_node = mesh->face_nodes[mesh->face_node_offsets[f]+n];
        int new_node = node_map[old_node];
        cropped_mesh->face_nodes[cropped_mesh->face_node_offsets[new_face]+n] = new_node;
      }
    }
  }

  // Copy node positions.
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    if (node_map[n] != -1)
      cropped_mesh->nodes[node_map[n]] = mesh->nodes[n];
  }

  // Build the rest.
  mesh_construct_edges(cropped_mesh);
  mesh_compute_geometry(cropped_mesh);
  
  // Create the boundary faces tag.
  int* bf_tag = mesh_create_tag(cropped_mesh->face_tags, sp_func_name(boundary_func), boundary_faces->size);
  int pos = 0, i = 0, face;
  while (int_unordered_set_next(boundary_faces, &pos, &face))
    bf_tag[i++] = face;

  // Clean up.
  polymec_free(node_map);
  polymec_free(face_map);
  int_unordered_set_free(boundary_faces);
  int_unordered_set_free(outside_cells);

  return cropped_mesh;
}

