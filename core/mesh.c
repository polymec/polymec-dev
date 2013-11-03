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

#include <stdlib.h>
#include "core/mesh.h"
#include "core/mesh_storage.h"
#include "core/unordered_set.h"
#include "core/linear_algebra.h"

// Generic tagging functions -- defined in tagger.c.
extern tagger_t* tagger_new(ARENA* arena);
extern void tagger_free(tagger_t* tagger);
extern int* tagger_create_tag(tagger_t* tagger, const char* tag_name, int size);
extern int* tagger_tag(tagger_t* tagger, const char* tag_name, int* size);
extern bool tagger_has_tag(tagger_t* tagger, const char* tag_name);
extern void tagger_delete_tag(tagger_t* tagger, const char* tag_name);
extern bool tagger_set_property(tagger_t* tagger, const char* tag_name, const char* property_name, void* data, void (*dtor)(void*));
extern void* tagger_property(tagger_t* tagger, const char* tag_name, const char* property_name);
extern void tagger_delete_property(tagger_t* tagger, const char* tag_name, const char* property_name);
extern void tagger_rename_tag(tagger_t* tagger, const char* old_tag_name, const char* new_tag_name);
extern bool tagger_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, int* tag_size);

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

mesh_t* mesh_new(int num_cells, int num_ghost_cells, int num_faces,
                 int num_edges, int num_nodes)
{
  ARENA* a = arena_open(&arena_defaults, 0);
  mesh_t* mesh = mesh_new_with_arena(a, num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);
  mesh->close_arena = true;
  return mesh;
}

mesh_t* mesh_new_with_arena(ARENA* arena, int num_cells, int num_ghost_cells, int num_faces,
                            int num_edges, int num_nodes)
{
  ASSERT(num_cells >= 0);
  ASSERT(num_ghost_cells >= 0);
  ASSERT(num_faces >= 0);
  ASSERT(num_edges >= 0);
  ASSERT(num_nodes >= 0);

  mesh_t* mesh = ARENA_MALLOC(arena, sizeof(mesh_t), 0);
  mesh->arena = arena;
  mesh->close_arena = false;

  // NOTE: We round stored elements up to the nearest power of 2.

  // Allocate cell information.
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;
  mesh->cell_face_offsets = ARENA_MALLOC(mesh->arena, sizeof(int)*(num_cells+num_ghost_cells+1), 0);
  memset(mesh->cell_face_offsets, 0, sizeof(int)*(num_cells+num_ghost_cells+1));
  int cell_face_cap = round_to_pow2(12 * (num_cells + num_ghost_cells));
  mesh->cell_faces = ARENA_MALLOC(mesh->arena, sizeof(int)*(cell_face_cap), 0);
  memset(mesh->cell_faces, 0, sizeof(int)*cell_face_cap);

  // Allocate face information.
  mesh->num_faces = num_faces;

  mesh->face_node_offsets = ARENA_MALLOC(mesh->arena, sizeof(int)*(num_faces+1), 0);
  memset(mesh->face_node_offsets, 0, sizeof(int)*(num_faces+1));
  int face_node_cap = round_to_pow2(6 * num_faces);
  mesh->face_nodes = ARENA_MALLOC(mesh->arena, sizeof(int)*(face_node_cap), 0);
  memset(mesh->face_nodes, 0, sizeof(int)*face_node_cap);

  mesh->face_edge_offsets = ARENA_MALLOC(mesh->arena, sizeof(int)*(num_faces+1), 0);
  memset(mesh->face_edge_offsets, 0, sizeof(int)*(num_faces+1));
  int face_edge_cap = round_to_pow2(18 * num_faces);
  mesh->face_edges = ARENA_MALLOC(mesh->arena, sizeof(int)*(face_edge_cap), 0);
  memset(mesh->face_edges, 0, sizeof(int)*face_edge_cap);

  mesh->face_cells = ARENA_MALLOC(mesh->arena, sizeof(int)*2*num_faces, 0);
  memset(mesh->face_cells, 0, sizeof(int)*2*num_faces);

  // Allocate edge information.
  mesh->num_edges = num_edges;
  mesh->edge_nodes = ARENA_MALLOC(mesh->arena, sizeof(int)*2*num_edges, 0);
  memset(mesh->edge_nodes, 0, sizeof(int)*2*num_edges);

  // Allocate node information.
  mesh->num_nodes = num_nodes;
  mesh->nodes = ARENA_MALLOC(mesh->arena, sizeof(point_t)*num_nodes, 0);
  memset(mesh->nodes, 0, sizeof(point_t)*num_nodes);

  // Allocate geometric data.
  mesh->cell_volumes = ARENA_MALLOC(mesh->arena, sizeof(double)*(num_cells+num_ghost_cells), 0);
  mesh->cell_centers = ARENA_MALLOC(mesh->arena, sizeof(point_t)*(num_cells+num_ghost_cells), 0);
  mesh->face_centers = ARENA_MALLOC(mesh->arena, sizeof(point_t)*num_faces, 0);
  mesh->face_areas = ARENA_MALLOC(mesh->arena, sizeof(double)*num_faces, 0);
  mesh->face_normals = ARENA_MALLOC(mesh->arena, sizeof(vector_t)*num_faces, 0);

  // Storage information.
  mesh->storage = mesh_storage_new_with_arena(arena);
  mesh->storage->cell_face_capacity = cell_face_cap;
  mesh->storage->face_node_capacity = face_node_cap;
  mesh->storage->face_edge_capacity = face_edge_cap;

  // Allocate tagging mechanisms.
  mesh->cell_tags = tagger_new(mesh->arena);
  mesh->face_tags = tagger_new(mesh->arena);
  mesh->edge_tags = tagger_new(mesh->arena);
  mesh->node_tags = tagger_new(mesh->arena);

  // Now we create a bogus tag that we can use to store mesh properties.
  int* prop_tag = mesh_create_tag(mesh->cell_tags, "properties", 1);
  prop_tag[0] = 0;

  return mesh;
}

void mesh_free(mesh_t* mesh)
{
  ASSERT(mesh != NULL);

  // Destroy tags.
  tagger_free(mesh->node_tags);
  tagger_free(mesh->edge_tags);
  tagger_free(mesh->face_tags);
  tagger_free(mesh->cell_tags);

  // Destroy storage metadata.
  mesh_storage_free(mesh->storage);

  // Destroy geometric information.
  ARENA_FREE(mesh->arena, mesh->face_normals);
  ARENA_FREE(mesh->arena, mesh->face_areas);
  ARENA_FREE(mesh->arena, mesh->face_centers);
  ARENA_FREE(mesh->arena, mesh->cell_centers);
  ARENA_FREE(mesh->arena, mesh->cell_volumes);

  // Destroy nodes.
  ARENA_FREE(mesh->arena, mesh->nodes);

  // Destroy connectivity.
  ARENA_FREE(mesh->arena, mesh->edge_nodes);
  ARENA_FREE(mesh->arena, mesh->face_cells);
  ARENA_FREE(mesh->arena, mesh->face_edges);
  ARENA_FREE(mesh->arena, mesh->face_edge_offsets);
  ARENA_FREE(mesh->arena, mesh->face_nodes);
  ARENA_FREE(mesh->arena, mesh->face_node_offsets);
  ARENA_FREE(mesh->arena, mesh->cell_faces);

  // Destroy the mesh itself (and possibly the arena).
  ARENA* arena = mesh->arena;
  bool close_arena = mesh->close_arena;
  ARENA_FREE(arena, mesh);
  if (close_arena)
    arena_close(arena);
}

void mesh_verify(mesh_t* mesh)
{
  // All cells must have at least 4 faces.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (mesh_cell_num_faces(mesh, c) < 4)
    {
      polymec_error("mesh_check_consistency: polyhedral cell %d has only %d faces.", 
                    c, mesh_cell_num_faces(mesh, c));
    }
  }

  // All faces must have at least 3 nodes/edges.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int ne = mesh_face_num_edges(mesh, f);
    if (ne == 0)
      polymec_error("mesh_check_consistency: polygonal face %d has no edges!", f);
    if (ne < 3)
      polymec_error("mesh_check_consistency: polygonal face %d has only %d edges.", f, ne);
  }

  // Check cell-face topology.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    int pos = 0, f;
    while (mesh_next_cell_face(mesh, c, &pos, &f))
    {
      if ((mesh->face_cells[2*f] != c) && (mesh->face_cells[2*f+1] != c))
        polymec_error("cell %d has face %d but is not attached to it.", c, f);
    }
  }
}

void mesh_set_property(mesh_t* mesh, const char* property, void* data, void (*dtor)(void*))
{
  // Use the bogus tag to store our junk.
  tagger_set_property(mesh->cell_tags, "properties", property, data, dtor);
}

void* mesh_property(mesh_t* mesh, const char* property)
{
  // Get this property from our bogus tag.
  return tagger_property(mesh->cell_tags, "properties", property);
}

void mesh_delete_property(mesh_t* mesh, const char* property)
{
  return tagger_delete_property(mesh->cell_tags, "properties", property);
}

int* mesh_create_tag(tagger_t* tagger, const char* tag, int num_indices)
{
  return tagger_create_tag(tagger, tag, num_indices);
}

int* mesh_tag(tagger_t* tagger, const char* tag, int* num_indices)
{
  return tagger_tag(tagger, tag, num_indices);
}

bool mesh_has_tag(tagger_t* tagger, const char* tag)
{
  return tagger_has_tag(tagger, tag);
}

bool mesh_tag_set_property(tagger_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*))
{
  return tagger_set_property(tagger, tag, property, data, destructor);
}

void* mesh_tag_property(tagger_t* tagger, const char* tag, const char* property)
{
  return tagger_property(tagger, tag, property);
}

void mesh_tag_delete_property(tagger_t* tagger, const char* tag, const char* property)
{
  tagger_delete_property(tagger, tag, property);
}

void mesh_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag)
{
  tagger_rename_tag(tagger, old_tag, new_tag);
}

void mesh_delete_tag(tagger_t* tagger, const char* tag)
{
  tagger_delete_tag(tagger, tag);
}

bool mesh_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, int* tag_size)
{
  return tagger_next_tag(tagger, pos, tag_name, tag_indices, tag_size);
}

void mesh_compute_geometry(mesh_t* mesh)
{
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    // Make sure each cell has at least 4 faces.
    ASSERT((mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell]) >= 4);

    // Compute cell centers and face centers for the cell, 
    // knowing that it's convex.
    mesh->cell_centers[cell].x = mesh->cell_centers[cell].y = mesh->cell_centers[cell].z = 0.0;
    int num_cell_nodes = 0, num_cell_faces = 0;
    int pos = 0, face;
    while (mesh_next_cell_face(mesh, cell, &pos, &face))
    {
      // NOTE: Only the primal cell of a face computes its center.
      if (cell == mesh->face_cells[2*face])
        mesh->face_centers[face].x = mesh->face_centers[face].y = mesh->face_centers[face].z = 0.0;
      int npos = 0, node;
      while (mesh_next_face_node(mesh, face, &npos, &node))
      {
        mesh->cell_centers[cell].x += mesh->nodes[node].x;
        mesh->cell_centers[cell].y += mesh->nodes[node].y;
        mesh->cell_centers[cell].z += mesh->nodes[node].z;
        if (cell == mesh->face_cells[2*face])
        {
          mesh->face_centers[face].x += mesh->nodes[node].x;
          mesh->face_centers[face].y += mesh->nodes[node].y;
          mesh->face_centers[face].z += mesh->nodes[node].z;
        }
      }
      if (cell == mesh->face_cells[2*face])
      {
        int nn = mesh_face_num_nodes(mesh, face);
        mesh->face_centers[face].x /= nn;
        mesh->face_centers[face].y /= nn;
        mesh->face_centers[face].z /= nn;
      }
      num_cell_nodes += mesh_face_num_nodes(mesh, face);
      ++num_cell_faces;
    }
    mesh->cell_centers[cell].x /= num_cell_nodes;
    mesh->cell_centers[cell].y /= num_cell_nodes;
    mesh->cell_centers[cell].z /= num_cell_nodes;

    // Use the preceding geometry to compute face areas and the 
    // cell's volume.
    mesh->cell_volumes[cell] = 0.0;
    pos = 0;
    while (mesh_next_cell_face(mesh, cell, &pos, &face))
    {
      double face_area = 0.0;
      vector_t face_normal = {0.0, 0.0, 0.0};
      int epos = 0, edge;
      while (mesh_next_face_edge(mesh, face, &epos, &edge))
      {
        // Construct a tetrahedron whose vertices are the cell center, 
        // the face center, and the two nodes of this edge. The volume 
        // of this tetrahedron contributes to the cell volume.
        vector_t v1, v2, v3, v2xv3;
        point_displacement(&mesh->face_centers[face], &mesh->cell_centers[cell], &v1);
        point_t xn1 = mesh->nodes[mesh->edge_nodes[2*edge]];
        point_t xn2 = mesh->nodes[mesh->edge_nodes[2*edge+1]];
        point_displacement(&mesh->face_centers[face], &xn1, &v2);
        point_displacement(&mesh->face_centers[face], &xn2, &v3);
        vector_cross(&v2, &v3, &v2xv3);
        double tet_volume = fabs(vector_dot(&v1, &v2xv3))/6.0;
        mesh->cell_volumes[cell] += tet_volume;

        // Now take the face of the tet whose vertices are the face center 
        // and the two nodes. The area of this tet contributes to the 
        // face's area.
        double tri_area = 0.5*vector_mag(&v2xv3);
        face_area += tri_area;
        face_normal = v2xv3;
      }
      // Only the primal cell of a face computes its area/normal.
      if (cell == mesh->face_cells[2*face])
      {
        vector_normalize(&face_normal);
        mesh->face_normals[face] = face_normal;
        mesh->face_areas[face] = face_area;
      }
    }
  }
}

