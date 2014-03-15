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

#include <stdlib.h>
#include "core/mesh.h"
#include "core/mesh_storage.h"
#include "core/unordered_set.h"
#include "core/table.h"

// Mesh features.
const char* TETRAHEDRAL = "tetrahedral";

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

mesh_t* mesh_new(MPI_Comm comm, int num_cells, int num_ghost_cells, 
                 int num_faces, int num_nodes)
{
  ARENA* a = arena_open(&arena_defaults, 0);
  mesh_t* mesh = mesh_new_with_arena(a, comm, num_cells, num_ghost_cells, num_faces, num_nodes);
  mesh->close_arena = true;
  return mesh;
}

mesh_t* mesh_new_with_arena(ARENA* arena, MPI_Comm comm, int num_cells, 
                            int num_ghost_cells, int num_faces, int num_nodes)
{
  ASSERT(num_cells >= 0);
  ASSERT(num_ghost_cells >= 0);
  ASSERT(num_faces >= 0);
  ASSERT(num_nodes >= 0);

  mesh_t* mesh = ARENA_MALLOC(arena, sizeof(mesh_t), 0);
  mesh->arena = arena;
  mesh->comm = comm;
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
  int face_edge_cap = round_to_pow2(4 * num_faces);
  mesh->face_edges = ARENA_MALLOC(mesh->arena, sizeof(int)*(face_edge_cap), 0);
  memset(mesh->face_edges, 0, sizeof(int)*face_edge_cap);

  mesh->face_cells = ARENA_MALLOC(mesh->arena, sizeof(int)*2*num_faces, 0);
  for (int f = 0; f < 2*mesh->num_faces; ++f)
    mesh->face_cells[f] = -1;

  // Allocate edge information.
  mesh->num_edges = 0;
  mesh->edge_nodes = NULL;

  // Allocate node information.
  mesh->num_nodes = num_nodes;
  mesh->nodes = ARENA_MALLOC(mesh->arena, sizeof(point_t)*num_nodes, 0);
  memset(mesh->nodes, 0, sizeof(point_t)*num_nodes);

  // Allocate geometric data.
  mesh->cell_volumes = ARENA_MALLOC(mesh->arena, sizeof(real_t)*(num_cells+num_ghost_cells), 0);
  mesh->cell_centers = ARENA_MALLOC(mesh->arena, sizeof(point_t)*(num_cells+num_ghost_cells), 0);
  mesh->face_centers = ARENA_MALLOC(mesh->arena, sizeof(point_t)*num_faces, 0);
  mesh->face_areas = ARENA_MALLOC(mesh->arena, sizeof(real_t)*num_faces, 0);
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

  // We also create a bogus tag that we can use to store mesh features.
  int* feature_tag = mesh_create_tag(mesh->cell_tags, "features", 1);
  feature_tag[0] = 0;

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
    while (mesh_cell_next_face(mesh, c, &pos, &f))
    {
      if ((mesh->face_cells[2*f] != c) && (mesh->face_cells[2*f+1] != c))
        polymec_error("cell %d has face %d but is not attached to it.", c, f);
    }
  }
}

mesh_t* mesh_clone(mesh_t* mesh)
{
  mesh_t* clone = mesh_new(MPI_COMM_WORLD, mesh->num_cells, mesh->num_ghost_cells,
                           mesh->num_faces, mesh->num_nodes);

  // Cell stuff.
  memcpy(clone->cell_face_offsets, mesh->cell_face_offsets, sizeof(int)*(mesh->num_cells+1));
  int num_cell_faces = clone->cell_face_offsets[clone->num_cells];
  if (clone->storage->cell_face_capacity < num_cell_faces)
  {
    clone->cell_faces = ARENA_REALLOC(clone->arena, clone->cell_faces, sizeof(int)*num_cell_faces, 0);
    clone->storage->cell_face_capacity = num_cell_faces;
  }
  memcpy(clone->cell_faces, mesh->cell_faces, sizeof(int)*num_cell_faces);

  // Face stuff.
  memcpy(clone->face_node_offsets, mesh->face_node_offsets, sizeof(int)*(mesh->num_faces+1));
  int num_face_nodes = clone->face_node_offsets[clone->num_faces];
  if (clone->storage->face_node_capacity < num_face_nodes)
  {
    clone->face_nodes = ARENA_REALLOC(clone->arena, clone->face_nodes, sizeof(int)*num_face_nodes, 0);
    clone->storage->face_node_capacity = num_face_nodes;
  }
  memcpy(clone->face_nodes, mesh->face_nodes, sizeof(int)*num_face_nodes);
  memcpy(clone->face_cells, mesh->face_cells, sizeof(int)*2*clone->num_faces);

  // Node stuff.
  memcpy(clone->nodes, mesh->nodes, sizeof(point_t)*clone->num_nodes);

  // Construct edges.
  mesh_construct_edges(clone);

  // Geometry stuff.
  memcpy(clone->cell_volumes, mesh->cell_volumes, sizeof(real_t)*clone->num_cells);
  memcpy(clone->cell_centers, mesh->cell_centers, sizeof(point_t)*clone->num_cells);
  memcpy(clone->face_centers, mesh->face_centers, sizeof(point_t)*clone->num_faces);
  memcpy(clone->face_areas, mesh->face_areas, sizeof(real_t)*clone->num_faces);
  memcpy(clone->face_normals, mesh->face_normals, sizeof(vector_t)*clone->num_faces);

  // Tags.
  tagger_copy(clone->node_tags, mesh->node_tags);
  tagger_copy(clone->edge_tags, mesh->edge_tags);
  tagger_copy(clone->face_tags, mesh->face_tags);
  tagger_copy(clone->cell_tags, mesh->cell_tags);

  return clone;
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

void mesh_add_feature(mesh_t* mesh, const char* feature)
{
  // Use the bogus tag to store our junk.
  bool* data = malloc(sizeof(bool));
  *data = true;
  tagger_set_property(mesh->cell_tags, "features", feature, data, free);
}

bool mesh_has_feature(mesh_t* mesh, const char* feature)
{
  // Get this feature from our bogus tag.
  return (tagger_property(mesh->cell_tags, "features", feature) != NULL);
}

void mesh_delete_feature(mesh_t* mesh, const char* feature)
{
  return tagger_delete_property(mesh->cell_tags, "features", feature);
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
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      ASSERT(face < mesh->num_faces);
      // NOTE: Only the primal cell of a face computes its center.
      if (cell == mesh->face_cells[2*face])
        mesh->face_centers[face].x = mesh->face_centers[face].y = mesh->face_centers[face].z = 0.0;
      int npos = 0, node;
      while (mesh_face_next_node(mesh, face, &npos, &node))
      {
        ASSERT(node >= 0);
        ASSERT(node < mesh->num_nodes);
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
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      real_t face_area = 0.0;
      vector_t face_normal = {0.0, 0.0, 0.0};
      int epos = 0, edge;
      while (mesh_face_next_edge(mesh, face, &epos, &edge))
      {
        ASSERT(edge >= 0);
        ASSERT(edge < mesh->num_edges);
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
        real_t tet_volume = fabs(vector_dot(&v1, &v2xv3))/6.0;
        mesh->cell_volumes[cell] += tet_volume;

        // Now take the face of the tet whose vertices are the face center 
        // and the two nodes. The area of this tet contributes to the 
        // face's area.
        real_t tri_area = 0.5*vector_mag(&v2xv3);
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

void mesh_construct_edges(mesh_t* mesh)
{
  ASSERT(mesh->num_edges == 0);
  ASSERT(mesh->edge_nodes == NULL);

  // Construct edge information.
  int_table_t* edge_for_nodes = int_table_new();
  {
    int num_edges = 0;
    for (int f = 0; f < mesh->num_faces; ++f)
    {
      int offset = 4*f;
      mesh->face_edge_offsets[f] = offset;
      for (int n = 0; n < 4; ++n)
      {
        // Make room.
        if (mesh->storage->face_edge_capacity <= offset+n)
        {
          while (mesh->storage->face_edge_capacity <= offset+n)
            mesh->storage->face_edge_capacity *= 2;
          mesh->face_edges = ARENA_REALLOC(mesh->arena, mesh->face_edges, sizeof(int) * mesh->storage->face_edge_capacity, 0);
        }

        int n1 = (int)mesh->face_nodes[offset+n];
        int n2 = (int)mesh->face_nodes[offset+(n+1)%4];
        if (!int_table_contains(edge_for_nodes, MIN(n1, n2), MAX(n1, n2)))
        {
          int_table_insert(edge_for_nodes, MIN(n1, n2), MAX(n1, n2), num_edges);
          mesh->face_edges[offset+n] = num_edges;
          ++num_edges;
        }
        else
          mesh->face_edges[offset+n] = *int_table_get(edge_for_nodes, MIN(n1, n2), MAX(n1, n2));
      }
    }
    if (num_edges > 0)
    {
      mesh->num_edges = num_edges;
      mesh->edge_nodes = ARENA_MALLOC(mesh->arena, 2 * sizeof(int) * mesh->num_edges, 0);
      int_table_cell_pos_t pos = int_table_start(edge_for_nodes);
      int n1, n2, e;
      while (int_table_next_cell(edge_for_nodes, &pos, &n1, &n2, &e))
      {
        mesh->edge_nodes[2*e] = n1;
        mesh->edge_nodes[2*e+1] = n2;
      }
    }
  }
  int_table_free(edge_for_nodes);
}

static size_t mesh_byte_size(void* obj)
{
  mesh_t* mesh = obj;
  
  size_t basic_storage = 
    // cell stuff
    2*sizeof(int) + (mesh->num_cells+1) * sizeof(int) + 
    mesh->cell_face_offsets[mesh->num_cells] * sizeof(int) + 
    // face stuff
    sizeof(int) + (mesh->num_faces+1) * sizeof(int) + 
    mesh->face_node_offsets[mesh->num_faces] * sizeof(int) + 
    2*mesh->num_faces * sizeof(int) + 
    // node stuff
    sizeof(int) + mesh->num_nodes*sizeof(point_t) + 
    // geometry
    mesh->num_cells*(sizeof(real_t) + sizeof(point_t)) + 
    mesh->num_faces*(sizeof(point_t) + sizeof(real_t) + sizeof(vector_t));
  
  size_t tag_storage = 4 * sizeof(int); // Numbers of tags.
  int pos = 0, *tag, tag_size;
  char* tag_name;
  while (mesh_next_tag(mesh->node_tags, &pos, &tag_name, &tag, &tag_size))
  {
    tag_storage += sizeof(int) + strlen(tag_name) * sizeof(char);
    tag_storage += sizeof(int) + tag_size * sizeof(int);
  }
  pos = 0;
  while (mesh_next_tag(mesh->edge_tags, &pos, &tag_name, &tag, &tag_size))
  {
    tag_storage += sizeof(int) + strlen(tag_name) * sizeof(char);
    tag_storage += sizeof(int) + tag_size * sizeof(int);
  }
  pos = 0;
  while (mesh_next_tag(mesh->face_tags, &pos, &tag_name, &tag, &tag_size))
  {
    tag_storage += sizeof(int) + strlen(tag_name) * sizeof(char);
    tag_storage += sizeof(int) + tag_size * sizeof(int);
  }
  pos = 0;
  while (mesh_next_tag(mesh->cell_tags, &pos, &tag_name, &tag, &tag_size))
  {
    tag_storage += sizeof(int) + strlen(tag_name) * sizeof(char);
    tag_storage += sizeof(int) + tag_size * sizeof(int);
  }

  return basic_storage + tag_storage;
}

static void byte_array_read_tags(byte_array_t* bytes, size_t* offset, tagger_t* tagger)
{
  int num_tags;
  byte_array_read_ints(bytes, 1, offset, &num_tags);
  for (int i = 0; i < num_tags; ++i)
  {
    int tag_name_len;
    byte_array_read_ints(bytes, 1, offset, &tag_name_len);
    char tag_name[tag_name_len+1];
    byte_array_read_chars(bytes, tag_name_len, offset, tag_name);
    tag_name[tag_name_len] = '\0';
    int tag_size;
    byte_array_read_ints(bytes, 1, offset, &tag_size);
    mesh_delete_tag(tagger, tag_name);
    int* tag = mesh_create_tag(tagger, tag_name, tag_size);
    byte_array_read_ints(bytes, tag_size, offset, tag);
  }
}

static void* mesh_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the number of cells, faces, nodes, and allocate a mesh
  // accordingly.
  int num_cells, num_ghost_cells, num_faces, num_nodes;
  byte_array_read_ints(bytes, 1, offset, &num_cells);
  byte_array_read_ints(bytes, 1, offset, &num_ghost_cells);
  byte_array_read_ints(bytes, 1, offset, &num_faces);
  byte_array_read_ints(bytes, 1, offset, &num_nodes);

  mesh_t* mesh = mesh_new(MPI_COMM_WORLD, num_cells, num_ghost_cells,
                          num_faces, num_nodes);

  // Read all the cell stuff.
  byte_array_read_ints(bytes, num_cells+1, offset, mesh->cell_face_offsets);
  int num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
  if (mesh->storage->cell_face_capacity < num_cell_faces)
  {
    mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int)*num_cell_faces, 0);
    mesh->storage->cell_face_capacity = num_cell_faces;
  }
  byte_array_read_ints(bytes, num_cell_faces, offset, mesh->cell_faces);

  // Face stuff.
  byte_array_read_ints(bytes, num_faces+1, offset, mesh->face_node_offsets);
  int num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
  if (mesh->storage->face_node_capacity < num_face_nodes)
  {
    mesh->face_nodes = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int)*num_face_nodes, 0);
    mesh->storage->face_node_capacity = num_face_nodes;
  }
  byte_array_read_ints(bytes, num_face_nodes, offset, mesh->face_nodes);
  byte_array_read_ints(bytes, 2*num_faces, offset, mesh->face_cells);

  // Node stuff.
  byte_array_read_points(bytes, num_nodes, offset, mesh->nodes);

  // Construct edges.
  mesh_construct_edges(mesh);

  // Geometry stuff.
  byte_array_read_reals(bytes, num_cells, offset, mesh->cell_volumes);
  byte_array_read_points(bytes, num_cells, offset, mesh->cell_centers);
  byte_array_read_points(bytes, num_faces, offset, mesh->face_centers);
  byte_array_read_reals(bytes, num_faces, offset, mesh->face_areas);
  byte_array_read_vectors(bytes, num_faces, offset, mesh->face_normals);

  // Tag stuff.
  byte_array_read_tags(bytes, offset, mesh->cell_tags);
  byte_array_read_tags(bytes, offset, mesh->face_tags);
  byte_array_read_tags(bytes, offset, mesh->edge_tags);
  byte_array_read_tags(bytes, offset, mesh->node_tags);

  return mesh;
}

static void byte_array_write_tags(byte_array_t* bytes, tagger_t* tagger, size_t* offset)
{
  // Count up the tags.
  int pos = 0, *tag, tag_size, num_tags = 0;
  char* tag_name;
  while (mesh_next_tag(tagger, &pos, &tag_name, &tag, &tag_size))
    ++num_tags;
  byte_array_write_ints(bytes, 1, &num_tags, offset);

  pos = 0;
  while (mesh_next_tag(tagger, &pos, &tag_name, &tag, &tag_size))
  {
    int tag_name_len = strlen(tag_name);
    byte_array_write_ints(bytes, 1, &tag_name_len, offset);
    byte_array_write_chars(bytes, tag_name_len, tag_name, offset);
    byte_array_write_ints(bytes, 1, &tag_size, offset);
    byte_array_write_ints(bytes, tag_size, tag, offset);
  }
}

static void mesh_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  mesh_t* mesh = obj;

  // Write the number of cells, faces, nodes.
  byte_array_write_ints(bytes, 1, &mesh->num_cells, offset);
  byte_array_write_ints(bytes, 1, &mesh->num_ghost_cells, offset);
  byte_array_write_ints(bytes, 1, &mesh->num_faces, offset);
  byte_array_write_ints(bytes, 1, &mesh->num_nodes, offset);

  // Write all the cell stuff.
  byte_array_write_ints(bytes, mesh->num_cells+1, mesh->cell_face_offsets, offset);
  byte_array_write_ints(bytes, mesh->cell_face_offsets[mesh->num_cells], mesh->cell_faces, offset);

  // Face stuff.
  byte_array_write_ints(bytes, mesh->num_faces+1, mesh->face_node_offsets, offset);
  byte_array_write_ints(bytes, mesh->face_node_offsets[mesh->num_faces], mesh->face_nodes, offset);

  byte_array_write_ints(bytes, 2*mesh->num_faces, mesh->face_cells, offset);

  // Node stuff.
  byte_array_write_points(bytes, mesh->num_nodes, mesh->nodes, offset);

  // Geometry stuff.
  byte_array_write_reals(bytes, mesh->num_cells, mesh->cell_volumes, offset);
  byte_array_write_points(bytes, mesh->num_cells, mesh->cell_centers, offset);
  byte_array_write_points(bytes, mesh->num_faces, mesh->face_centers, offset);
  byte_array_write_reals(bytes, mesh->num_faces, mesh->face_areas, offset);
  byte_array_write_vectors(bytes, mesh->num_faces, mesh->face_normals, offset);

  // Tag stuff.
  byte_array_write_tags(bytes, mesh->cell_tags, offset);
  byte_array_write_tags(bytes, mesh->face_tags, offset);
  byte_array_write_tags(bytes, mesh->edge_tags, offset);
  byte_array_write_tags(bytes, mesh->node_tags, offset);
}

serializer_t* mesh_serializer()
{
  return serializer_new(mesh_byte_size, mesh_byte_read, mesh_byte_write);
}

adj_graph_t* graph_from_mesh_cells(mesh_t* mesh)
{
  // Create a graph whose vertices are the mesh's cells.
  adj_graph_t* g = adj_graph_new(mesh->comm, mesh->num_cells);

  // Allocate space in the graph for the edges (faces connecting cells).
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    // How many faces don't have opposite cells?
    int outer_faces = 0;
    for (int j = mesh->cell_face_offsets[i]; j < mesh->cell_face_offsets[i+1]; ++j)
    {
      int f = mesh->cell_faces[j];
      if (f < 0) f = ~f;
      if (mesh->face_cells[2*f+1] == -1)
        ++outer_faces;
    }
    adj_graph_set_num_edges(g, i, mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i] - outer_faces);
  }

  // Now fill in the edges.
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int* edges = adj_graph_edges(g, i);
    int offset = 0;
    for (int j = mesh->cell_face_offsets[i]; j < mesh->cell_face_offsets[i+1]; ++j)
    {
      int f = mesh->cell_faces[j];
      if (f < 0) f = ~f;
      if (mesh->face_cells[2*f+1] != -1)
      {
        int c = (i == mesh->face_cells[2*f]) ? mesh->face_cells[2*f+1] : mesh->face_cells[2*f];
        edges[offset] = c;
        ++offset;
      }
    }
  }

  return g;
}

