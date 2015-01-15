// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/mesh.h"
#include "core/table.h"
#include "core/unordered_set.h"

// Mesh features.
const char* TETRAHEDRAL = "tetrahedral";

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

// Mesh storage stuff.
struct mesh_storage_t
{
  int cell_face_capacity;
  int face_edge_capacity;
  int face_node_capacity;
  exchanger_t* exchanger;
};

// Initializes a new storage mechanism for a mesh.
static mesh_storage_t* mesh_storage_new(MPI_Comm comm)
{
  mesh_storage_t* storage = polymec_malloc(sizeof(mesh_storage_t));
  storage->cell_face_capacity = 0;
  storage->face_edge_capacity = 0;
  storage->face_node_capacity = 0;
  storage->exchanger = exchanger_new(comm);
  return storage;
}

// Frees the given storage mechanism.
static void mesh_storage_free(mesh_storage_t* storage)
{
  exchanger_free(storage->exchanger);
  polymec_free(storage);
}

mesh_t* mesh_new(MPI_Comm comm, int num_cells, int num_ghost_cells, 
                 int num_faces, int num_nodes)
{
  ASSERT(num_cells >= 0);
  ASSERT(num_ghost_cells >= 0);
  ASSERT(num_faces >= 0);
  ASSERT(num_nodes >= 0);

  mesh_t* mesh = polymec_malloc(sizeof(mesh_t));
  mesh->comm = comm;

  // NOTE: We round stored elements up to the nearest power of 2.

  // Allocate cell information.
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;
  mesh->cell_face_offsets = polymec_malloc(sizeof(int)*(num_cells+1));
  memset(mesh->cell_face_offsets, 0, sizeof(int)*(num_cells+1));
  int cell_face_cap = round_to_pow2(12 * num_cells);
  mesh->cell_faces = polymec_malloc(sizeof(int)*(cell_face_cap));
  memset(mesh->cell_faces, 0, sizeof(int)*cell_face_cap);

  // Allocate face information.
  mesh->num_faces = num_faces;

  mesh->face_node_offsets = polymec_malloc(sizeof(int)*(num_faces+1));
  memset(mesh->face_node_offsets, 0, sizeof(int)*(num_faces+1));
  int face_node_cap = round_to_pow2(6 * num_faces);
  mesh->face_nodes = polymec_malloc(sizeof(int)*(face_node_cap));
  memset(mesh->face_nodes, 0, sizeof(int)*face_node_cap);

  mesh->face_edge_offsets = polymec_malloc(sizeof(int)*(num_faces+1));
  memset(mesh->face_edge_offsets, 0, sizeof(int)*(num_faces+1));
  int face_edge_cap = round_to_pow2(4 * num_faces);
  mesh->face_edges = polymec_malloc(sizeof(int)*(face_edge_cap));
  memset(mesh->face_edges, 0, sizeof(int)*face_edge_cap);

  mesh->face_cells = polymec_malloc(sizeof(int)*2*num_faces);
  for (int f = 0; f < 2*mesh->num_faces; ++f)
    mesh->face_cells[f] = -1;

  // Allocate edge information.
  mesh->num_edges = 0;
  mesh->edge_nodes = NULL;

  // Allocate node information.
  mesh->num_nodes = num_nodes;
  mesh->nodes = polymec_malloc(sizeof(point_t)*num_nodes);
  memset(mesh->nodes, 0, sizeof(point_t)*num_nodes);

  // Allocate geometric data.
  int total_num_cells = num_cells + num_ghost_cells;
  mesh->cell_volumes = polymec_malloc(sizeof(real_t)*total_num_cells);
  mesh->cell_centers = polymec_malloc(sizeof(point_t)*total_num_cells);
  mesh->face_centers = polymec_malloc(sizeof(point_t)*num_faces);
  mesh->face_areas = polymec_malloc(sizeof(real_t)*num_faces);
  mesh->face_normals = polymec_malloc(sizeof(vector_t)*num_faces);

  // Storage information.
  mesh->storage = mesh_storage_new(mesh->comm);
  mesh->storage->cell_face_capacity = cell_face_cap;
  mesh->storage->face_node_capacity = face_node_cap;
  mesh->storage->face_edge_capacity = face_edge_cap;

  // Allocate tagging mechanisms.
  mesh->cell_tags = tagger_new();
  mesh->face_tags = tagger_new();
  mesh->edge_tags = tagger_new();
  mesh->node_tags = tagger_new();

  // Now we create a bogus tag that we can use to store mesh properties.
  int* prop_tag = mesh_create_tag(mesh->cell_tags, "properties", 1);
  prop_tag[0] = 0;

  // We also create a bogus tag that we can use to store mesh features.
  int* feature_tag = mesh_create_tag(mesh->cell_tags, "features", 1);
  feature_tag[0] = 0;

  return mesh;
}

mesh_t* mesh_new_with_cell_type(MPI_Comm comm, int num_cells, 
                                int num_ghost_cells, int num_faces, 
                                int num_nodes, int num_faces_per_cell,
                                int num_nodes_per_face)
{
  mesh_t* mesh = mesh_new(comm, num_cells, num_ghost_cells, num_faces, num_nodes);

  // Set up connectivity metadata.
  for (int c = 0; c < mesh->num_cells+1; ++c)
    mesh->cell_face_offsets[c] = num_faces_per_cell*c;
  for (int f = 0; f < mesh->num_faces+1; ++f)
    mesh->face_node_offsets[f] = num_nodes_per_face*f;
  mesh_reserve_connectivity_storage(mesh);
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
  polymec_free(mesh->face_normals);
  polymec_free(mesh->face_areas);
  polymec_free(mesh->face_centers);
  polymec_free(mesh->cell_centers);
  polymec_free(mesh->cell_volumes);

  // Destroy nodes.
  polymec_free(mesh->nodes);

  // Destroy connectivity.
  polymec_free(mesh->edge_nodes);
  polymec_free(mesh->face_cells);
  polymec_free(mesh->face_edges);
  polymec_free(mesh->face_edge_offsets);
  polymec_free(mesh->face_nodes);
  polymec_free(mesh->face_node_offsets);
  polymec_free(mesh->cell_faces);
  polymec_free(mesh->cell_face_offsets);

  // Destroy the mesh itself.
  polymec_free(mesh);
}

bool mesh_verify_topology(mesh_t* mesh, void (*handler)(const char* format, ...))
{
  // All cells must have at least 4 faces.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (mesh_cell_num_faces(mesh, c) < 4)
    {
      handler("mesh_verify_topology: polyhedral cell %d has only %d faces.", 
              c, mesh_cell_num_faces(mesh, c));
      return false;
    }
  }

  // All faces must have at least 3 nodes/edges.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int ne = mesh_face_num_edges(mesh, f);
    if (ne == 0)
    {
      handler("mesh_verify_topology: polygonal face %d has no edges!", f);
      return false;
    }
    if (ne < 3)
    {
      handler("mesh_verify_topology: polygonal face %d has only %d edges.", f, ne);
      return false;
    }
  }

  // Check cell-face topology.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    int pos = 0, f;
    while (mesh_cell_next_face(mesh, c, &pos, &f))
    {
      if ((mesh->face_cells[2*f] != c) && (mesh->face_cells[2*f+1] != c))
      {
        handler("mesh_verify_topology: cell %d has face %d but is not attached to it.", c, f);
        return false;
      }
    }
  }
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int pos = 0, ff;
    bool found_face = false;
    while (mesh_cell_next_face(mesh, mesh->face_cells[2*f], &pos, &ff))
    {
      if (ff == f) 
      {
        found_face = true;
        break;
      }
    }
    if (!found_face)
    {
      handler("mesh_verify_topology: face %d has cell %d but is not attached to it.", f, mesh->face_cells[2*f]);
      return false;
    }
    if (mesh->face_cells[2*f+1] != -1)
    {
      while (mesh_cell_next_face(mesh, mesh->face_cells[2*f], &pos, &ff))
      {
        if (ff == f) 
        {
          found_face = true;
          break;
        }
      }
      if (!found_face)
      {
        handler("mesh_verify_topology: face %d has cell %d but is not attached to it.", f, mesh->face_cells[2*f+1]);
        return false;
      }
    }
  }

  return true;
}

mesh_t* mesh_clone(mesh_t* mesh)
{
  mesh_t* clone = mesh_new(MPI_COMM_WORLD, mesh->num_cells, mesh->num_ghost_cells,
                           mesh->num_faces, mesh->num_nodes);

  // Connectivity metadata.
  memcpy(clone->cell_face_offsets, mesh->cell_face_offsets, sizeof(int)*(mesh->num_cells+1));
  memcpy(clone->face_node_offsets, mesh->face_node_offsets, sizeof(int)*(mesh->num_faces+1));
  mesh_reserve_connectivity_storage(clone);

  // Actual connectivity.
  int num_cell_faces = clone->cell_face_offsets[clone->num_cells];
  memcpy(clone->cell_faces, mesh->cell_faces, sizeof(int)*num_cell_faces);
  int num_face_nodes = clone->face_node_offsets[clone->num_faces];
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

void mesh_set_property(mesh_t* mesh, 
                       const char* property, 
                       void* data, 
                       serializer_t* serializer)
{
  // Use the bogus tag to store our junk.
  tagger_set_property(mesh->cell_tags, "properties", property, data, serializer);
}

void* mesh_property(mesh_t* mesh, const char* property)
{
  // Get this property from our bogus tag.
  return tagger_property(mesh->cell_tags, "properties", property);
}

void mesh_delete_property(mesh_t* mesh, const char* property)
{
  tagger_delete_property(mesh->cell_tags, "properties", property);
}

bool mesh_next_property(mesh_t* mesh, int* pos, 
                        char** prop_name, void** prop_data, 
                        serializer_t** prop_serializer)
{
  return tagger_next_property(mesh->cell_tags, "properties", pos, prop_name, 
                              prop_data, prop_serializer);
}

exchanger_t* mesh_exchanger(mesh_t* mesh)
{
  return mesh->storage->exchanger;
}

void mesh_set_exchanger(mesh_t* mesh, exchanger_t* ex)
{
  ASSERT(ex != NULL);
  if (mesh->storage->exchanger != NULL)
    exchanger_free(mesh->storage->exchanger);
  mesh->storage->exchanger = ex;
}

void mesh_add_feature(mesh_t* mesh, const char* feature)
{
  // Use the bogus tag to store our junk.
  bool* data = polymec_malloc(sizeof(bool));
  *data = true;
  serializer_t* ser = string_serializer();
  tagger_set_property(mesh->cell_tags, "features", feature, data, ser);
}

bool mesh_has_feature(mesh_t* mesh, const char* feature)
{
  // Get this feature from our bogus tag.
  return (tagger_property(mesh->cell_tags, "features", feature) != NULL);
}

void mesh_delete_feature(mesh_t* mesh, const char* feature)
{
  tagger_delete_property(mesh->cell_tags, "features", feature);
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

bool mesh_tag_set_property(tagger_t* tagger, const char* tag, const char* property, void* data, serializer_t* serializer)
{
  return tagger_set_property(tagger, tag, property, data, serializer);
}

void* mesh_tag_property(tagger_t* tagger, const char* tag, const char* property)
{
  return tagger_property(tagger, tag, property);
}

void mesh_tag_delete_property(tagger_t* tagger, const char* tag, const char* property)
{
  tagger_delete_property(tagger, tag, property);
}

bool mesh_tag_next_property(tagger_t* tagger, const char* tag, int* pos, 
                            char** prop_name, void** prop_data, 
                            serializer_t** prop_serializer)
{
  return tagger_next_property(tagger, tag, pos, prop_name, prop_data,
                              prop_serializer);
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
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    int num_cell_nodes = 0, num_cell_faces = 0;
    int pos = 0, face;
    while (mesh_cell_next_oriented_face(mesh, cell, &pos, &face))
    {
      int actual_face = (face >= 0) ? face : ~face;
      ASSERT(actual_face < mesh->num_faces);
      point_t xf = {.x = 0.0, .y = 0.0, .z = 0.0};
      int npos = 0, node;
      while (mesh_face_next_node(mesh, actual_face, &npos, &node))
      {
        ASSERT(node >= 0);
        ASSERT(node < mesh->num_nodes);
        point_t* xn = &mesh->nodes[node];
        xc.x += xn->x;
        xc.y += xn->y;
        xc.z += xn->z;
        // NOTE: Only the primal cell of a face computes its center.
        if (cell == mesh->face_cells[2*actual_face])
        {
          xf.x += xn->x;
          xf.y += xn->y;
          xf.z += xn->z;
        }
      }
      if (cell == mesh->face_cells[2*actual_face])
      {
        int nn = mesh_face_num_nodes(mesh, actual_face);
        xf.x /= nn;
        xf.y /= nn;
        xf.z /= nn;
        mesh->face_centers[actual_face] = xf;
      }
      num_cell_nodes += mesh_face_num_nodes(mesh, actual_face);
      ++num_cell_faces;
    }
    xc.x /= num_cell_nodes;
    xc.y /= num_cell_nodes;
    xc.z /= num_cell_nodes;
    mesh->cell_centers[cell] = xc;

    // Use the preceding geometry to compute face areas and the 
    // cell's volume.
    mesh->cell_volumes[cell] = 0.0;
    pos = 0;
    while (mesh_cell_next_oriented_face(mesh, cell, &pos, &face))
    {
      int actual_face = (face >= 0) ? face : ~face;
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
        point_t* xf = &mesh->face_centers[actual_face];
        point_displacement(xf, &mesh->cell_centers[cell], &v1);
        point_t xn1 = mesh->nodes[mesh->edge_nodes[2*edge]];
        point_t xn2 = mesh->nodes[mesh->edge_nodes[2*edge+1]];
        point_displacement(xf, &xn1, &v2);
        point_displacement(xf, &xn2, &v3);
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

      // The cell that stores the face gets 
      // to take responsibility for the face normals/areas.
      if (cell == mesh->face_cells[2*actual_face])
      {
        mesh->face_areas[actual_face] = face_area;
        vector_normalize(&face_normal);

        // Flip the normal vector if we need to.
        // FIXME: We should revisit the above code so this isn't needed.
        vector_t outward;
        point_t* xf = &mesh->face_centers[actual_face];
        point_displacement(&mesh->cell_centers[cell], xf, &outward);
        real_t n_o_cf = vector_dot(&face_normal, &outward);
        if (n_o_cf < 0.0) 
        {
          if (face == actual_face)
            vector_scale(&face_normal, -1.0);
        }
        else if (n_o_cf > 0.0)
        {
          if (face != actual_face)
            vector_scale(&face_normal, -1.0);
        }
        mesh->face_normals[actual_face] = face_normal;
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
  memcpy(mesh->face_edge_offsets, mesh->face_node_offsets, sizeof(int) * (mesh->num_faces + 1));
  {
    int num_edges = 0;
    for (int f = 0; f < mesh->num_faces; ++f)
    {
      int offset = mesh->face_edge_offsets[f];
      int num_face_edges = mesh->face_edge_offsets[f+1] - offset;
      for (int n = 0; n < num_face_edges; ++n)
      {
        // Make room.
        if (mesh->storage->face_edge_capacity <= offset+n)
        {
          while (mesh->storage->face_edge_capacity <= offset+n)
            mesh->storage->face_edge_capacity *= 2;
          mesh->face_edges = polymec_realloc(mesh->face_edges, sizeof(int) * mesh->storage->face_edge_capacity);
        }

        int n1 = (int)mesh->face_nodes[offset+n];
        int n2 = (int)mesh->face_nodes[offset+(n+1)%num_face_edges];
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
      mesh->edge_nodes = polymec_malloc(2 * sizeof(int) * mesh->num_edges);
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

void mesh_reserve_connectivity_storage(mesh_t* mesh)
{
  // Make sure metadata is in order.
  int num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
  ASSERT(num_cell_faces >= 4*mesh->num_cells); 
  int num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
  ASSERT(num_face_nodes >= 3*mesh->num_faces); 

  if (mesh->storage->cell_face_capacity <= num_cell_faces)
  {
    while (mesh->storage->cell_face_capacity <= num_cell_faces)
      mesh->storage->cell_face_capacity *= 2;
    mesh->cell_faces = polymec_realloc(mesh->cell_faces, sizeof(int) * mesh->storage->cell_face_capacity);
  }
  if (mesh->storage->face_node_capacity <= num_face_nodes)
  {
    while (mesh->storage->face_node_capacity <= num_face_nodes)
      mesh->storage->face_node_capacity *= 2;
    mesh->face_nodes = polymec_realloc(mesh->face_nodes, sizeof(int) * mesh->storage->face_node_capacity);
  }
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
    sizeof(int) + mesh->num_nodes*sizeof(point_t);
  
  // Tag-related storage.
  serializer_t* tag_s = tagger_serializer();
  size_t tag_storage = serializer_size(tag_s, mesh->cell_tags) + 
                       serializer_size(tag_s, mesh->face_tags) + 
                       serializer_size(tag_s, mesh->edge_tags) + 
                       serializer_size(tag_s, mesh->node_tags);
  tag_s = NULL;

  // Exchanger-related storage.
  serializer_t* ex_s = exchanger_serializer();
  size_t ex_storage = serializer_size(ex_s, mesh_exchanger(mesh));
  ex_s = NULL;

  return basic_storage + tag_storage + ex_storage;
}

static void* mesh_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the number of cells, faces, nodes, and allocate a mesh
  // accordingly.
  int num_cells, num_ghost_cells, num_faces, num_nodes;
  byte_array_read_ints(bytes, 1, &num_cells, offset);
  byte_array_read_ints(bytes, 1, &num_ghost_cells, offset);
  byte_array_read_ints(bytes, 1, &num_faces, offset);
  byte_array_read_ints(bytes, 1, &num_nodes, offset);

  mesh_t* mesh = mesh_new(MPI_COMM_WORLD, num_cells, num_ghost_cells,
                          num_faces, num_nodes);

  // Read all the connectivity metadata.
  byte_array_read_ints(bytes, num_cells+1, mesh->cell_face_offsets, offset);
  byte_array_read_ints(bytes, num_faces+1, mesh->face_node_offsets, offset);

  // Make sure that connectivity storage is sufficient.
  mesh_reserve_connectivity_storage(mesh);

  // Actual connectivity data.
  int num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
  byte_array_read_ints(bytes, num_cell_faces, mesh->cell_faces, offset);
  int num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
  byte_array_read_ints(bytes, num_face_nodes, mesh->face_nodes, offset);
  byte_array_read_ints(bytes, 2*num_faces, mesh->face_cells, offset);

  // Node stuff.
  byte_array_read_points(bytes, num_nodes, mesh->nodes, offset);

  // Construct edges and compute geometry.
  mesh_construct_edges(mesh);
  mesh_compute_geometry(mesh);

  // Tag stuff.
  tagger_free(mesh->cell_tags);
  tagger_free(mesh->face_tags);
  tagger_free(mesh->edge_tags);
  tagger_free(mesh->node_tags);
  serializer_t* ser = tagger_serializer();
  mesh->cell_tags = serializer_read(ser, bytes, offset);
  mesh->face_tags = serializer_read(ser, bytes, offset);
  mesh->edge_tags = serializer_read(ser, bytes, offset);
  mesh->node_tags = serializer_read(ser, bytes, offset);

  // Storage/exchanger stuff.
  mesh_storage_t* storage = mesh->storage;
  byte_array_read_ints(bytes, 1, &storage->cell_face_capacity, offset);
  byte_array_read_ints(bytes, 1, &storage->face_edge_capacity, offset);
  byte_array_read_ints(bytes, 1, &storage->face_node_capacity, offset);
  exchanger_free(storage->exchanger);
  ser = exchanger_serializer();
  storage->exchanger = serializer_read(ser, bytes, offset);

  return mesh;
}

static void mesh_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  mesh_t* mesh = obj;

  // Write the number of cells, faces, nodes.
  byte_array_write_ints(bytes, 1, &mesh->num_cells, offset);
  byte_array_write_ints(bytes, 1, &mesh->num_ghost_cells, offset);
  byte_array_write_ints(bytes, 1, &mesh->num_faces, offset);
  byte_array_write_ints(bytes, 1, &mesh->num_nodes, offset);

  // Write all the connectivity metadata.
  byte_array_write_ints(bytes, mesh->num_cells+1, mesh->cell_face_offsets, offset);
  byte_array_write_ints(bytes, mesh->num_faces+1, mesh->face_node_offsets, offset);

  // Write all the actual connectivity data.
  byte_array_write_ints(bytes, mesh->cell_face_offsets[mesh->num_cells], mesh->cell_faces, offset);
  byte_array_write_ints(bytes, mesh->face_node_offsets[mesh->num_faces], mesh->face_nodes, offset);
  byte_array_write_ints(bytes, 2*mesh->num_faces, mesh->face_cells, offset);

  // Node stuff.
  byte_array_write_points(bytes, mesh->num_nodes, mesh->nodes, offset);

  // Tag stuff.
  serializer_t* ser = tagger_serializer();
  serializer_write(ser, mesh->cell_tags, bytes, offset);
  serializer_write(ser, mesh->face_tags, bytes, offset);
  serializer_write(ser, mesh->edge_tags, bytes, offset);
  serializer_write(ser, mesh->node_tags, bytes, offset);

  // Storage/exchanger stuff.
  mesh_storage_t* storage = mesh->storage;
  byte_array_write_ints(bytes, 1, &storage->cell_face_capacity, offset);
  byte_array_write_ints(bytes, 1, &storage->face_edge_capacity, offset);
  byte_array_write_ints(bytes, 1, &storage->face_node_capacity, offset);
  ser = exchanger_serializer();
  serializer_write(ser, mesh_exchanger(mesh), bytes, offset);
  ser = NULL;
}

serializer_t* mesh_serializer()
{
  return serializer_new("mesh", mesh_byte_size, mesh_byte_read, mesh_byte_write, NULL);
}

exchanger_t* mesh_face_exchanger_new(mesh_t* mesh)
{
  exchanger_t* ex = exchanger_new(mesh->comm);
#if POLYMEC_HAVE_MPI
  // Get the mesh cell exchanger.
  exchanger_t* cell_ex = mesh_exchanger(mesh);

  // Traverse the send cells and create send "faces."
  int_array_t* send_faces = int_array_new();
  int_array_t* receive_faces = int_array_new();
  int pos = 0, proc, *s_indices, s_size; 
  while (exchanger_next_send(cell_ex, &pos, &proc, &s_indices, &s_size))
  {
    int_array_clear(send_faces);
    int_array_clear(receive_faces);

    // Get the receive transaction corresponding to this send one.
    int *r_indices, r_size; 
    bool have_it = exchanger_get_receive(cell_ex, proc, &r_indices, &r_size);
    ASSERT(have_it);

    // Now find the face separating each pair of cells.
    ASSERT(r_size == s_size);
    for (int i = 0; i < s_size; ++i)
    {
      int s_cell = s_indices[i];
      int r_cell = r_indices[i];

      int fpos = 0, face;
      while (mesh_cell_next_face(mesh, s_cell, &fpos, &face))
      {
        int neighbor = mesh_face_opp_cell(mesh, face, s_cell);
        if (neighbor == r_cell)
        {
          int_array_append(send_faces, 2*face);
          int_array_append(receive_faces, 2*face+1);
          break;
        }
      }
    }
    ASSERT(send_faces->size == s_size);
    ASSERT(receive_faces->size == r_size);

    // Set up the exchange.
    if (send_faces->size > 0)
      exchanger_set_send(ex, proc, send_faces->data, send_faces->size, true);
    if (receive_faces->size > 0)
      exchanger_set_receive(ex, proc, receive_faces->data, receive_faces->size, true);
  }
  int_array_free(send_faces);
  int_array_free(receive_faces);
#endif
  return ex;
}

exchanger_t* mesh_node_exchanger_new(mesh_t* mesh)
{
  exchanger_t* ex = exchanger_new(mesh->comm);
#if POLYMEC_HAVE_MPI
  // Create a face exchanger.
  exchanger_t* face_ex = mesh_face_exchanger_new(mesh);

  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  //------------------------------------------------------------------------
  // We define the notion of a process "owning" a node thus: a process owns 
  // a node if it has the lowest rank of all processes on which the node 
  // is represented. The tricky thing about this definition is that a node 
  // represented on a process p can be owned by a process, q (say), that p 
  // does not interact with directly in the sense of face exchanges. This 
  // means we have to search for the owner in two sweeps.
  //------------------------------------------------------------------------

  // Traverse the send faces and make a list of nodes associated with 
  // these faces/processes.
  int num_neighbors = exchanger_num_sends(face_ex);
  int neighbors[num_neighbors];
  int_array_t* face_nodes[num_neighbors];
  int pos = 0, proc, *indices, size, p = 0;
  while (exchanger_next_send(face_ex, &pos, &proc, &indices, &size))
  {
    neighbors[p] = proc;
    int_unordered_set_t* face_node_set = int_unordered_set_new();
    for (int f = 0; f < size; ++f)
    {
      int face = indices[f]/2; // see docs on face exchanger!
      int npos = 0, node;
      while (mesh_face_next_node(mesh, face, &npos, &node))
        int_unordered_set_insert(face_node_set, node);
    }
    face_nodes[p] = int_array_new();
    int npos = 0, node;
    while (int_unordered_set_next(face_node_set, &npos, &node))
      int_array_append(face_nodes[p], node);
    int_unordered_set_free(face_node_set);
    ++p;
  }

  //------------------------------------------------------------------------
  // Create a consistent indexing of our nodes with each of our neighbors
  // by matching up pairs of locate/remote nodes.
  //------------------------------------------------------------------------

  MPI_Request requests[2*num_neighbors];
  MPI_Status statuses[2*num_neighbors];

  // First, gather the local and remote node positions.
  pos = 0, p = 0;
  point_t* remote_node_positions[num_neighbors];
  while (exchanger_next_receive(face_ex, &pos, &proc, &indices, &size))
  {
    int num_nodes = face_nodes[p]->size;
    remote_node_positions[p] = polymec_malloc(sizeof(point_t) * num_nodes);
    MPI_Irecv(remote_node_positions[p], 3*num_nodes, MPI_REAL_T, proc, 0, 
              mesh->comm, &requests[p]);
    ++p;
  }
  pos = 0, p = 0;
  point_t* local_node_positions[num_neighbors];
  while (exchanger_next_send(face_ex, &pos, &proc, &indices, &size))
  {
    int num_nodes = face_nodes[p]->size;
    local_node_positions[p] = polymec_malloc(sizeof(point_t) * num_nodes);
    for (int n = 0; n < num_nodes; ++n)
      local_node_positions[p][n] = mesh->nodes[face_nodes[p]->data[n]];
    MPI_Isend(local_node_positions[p], 3*num_nodes, MPI_REAL_T, proc, 0, 
              mesh->comm, &requests[p + num_neighbors]);
    ++p;
  }
  MPI_Waitall(2 * num_neighbors, requests, statuses);

  // Now match each remote node with a local node by picking the remote 
  // node that is closest. The ordering of the pairs will adhere to that 
  // of the process with the lower rank.
  pos = 0, p = 0;
  int_array_t* matched_nodes[num_neighbors];
  int_unordered_set_t* already_matched = int_unordered_set_new();
  while (exchanger_next_receive(face_ex, &pos, &proc, &indices, &size))
  {
    int num_nodes = face_nodes[p]->size;
    int_unordered_set_clear(already_matched);
    matched_nodes[p] = int_array_new();
    for (int n1 = 0; n1 < num_nodes; ++n1)
    {
      point_t* x1 = &(local_node_positions[p][n1]);
      real_t d_min = FLT_MAX;
      int i_min = -1;
      for (int n2 = 0; n2 < num_nodes; ++n2)
      {
        point_t* x2 = &(remote_node_positions[p][n2]);
        real_t d12 = point_distance(x1, x2);
        if (d12 < d_min)
        {
          d_min = d12;
          i_min = (rank < proc) ? n1 : n2;
        }
      }
      if (int_unordered_set_contains(already_matched, i_min))
      {
        polymec_error("mesh_node_exchanger_new: inconsistent node geometry!\n"
                      "  Cannot construct node exchanger.");
      }
      else
        int_unordered_set_insert(already_matched, i_min);
      int_array_append(matched_nodes[p], i_min);
    }
    ++p;
  }
  int_unordered_set_free(already_matched);

  // Preliminary cleanup.
  for (int p = 0; p < num_neighbors; ++p)
  {
    polymec_free(local_node_positions[p]);
    polymec_free(remote_node_positions[p]);
  }

  //------------------------------------------------------------------------
  // At this point, matched_nodes[p] contains an ordered array of node 
  // indices that can be used to construct send/receive arrays for this 
  // process and its pth neighbor. It remains to cull nodes from these 
  // arrays that are owned by a third process, and to construct separate 
  // arrays to express the connection between that third process and this one.
  //------------------------------------------------------------------------

  // We begin by figuring out which process owns each of the nodes we are 
  // "discussing."
  int_int_unordered_map_t* node_owners = int_int_unordered_map_new();
  for (int p = 0; p < num_neighbors; ++p)
  {
    int neighbor = neighbors[p];
    for (int n = 0; n < matched_nodes[p]->size; ++n)
    {
      int node = matched_nodes[p]->data[n];
      int owner = (neighbor < rank) ? neighbor : rank;
      int* owner_p = int_int_unordered_map_get(node_owners, node);
      if (owner_p == NULL)
        int_int_unordered_map_insert(node_owners, node, owner);
      else
      {
        if (owner < *owner_p)
          int_int_unordered_map_insert(node_owners, node, owner);
      }
    }
  }

  // Construct an owner array for each of the nodes.
  int_array_t* matched_node_owners[num_neighbors];
  for (int p = 0; p < num_neighbors; ++p)
  {
    matched_node_owners[p] = int_array_new();
    for (int n = 0; n < matched_nodes[p]->size; ++n)
    {
      int node = matched_nodes[p]->data[n];
      int owner = *int_int_unordered_map_get(node_owners, node);
      int_array_append(matched_node_owners[p], owner);
    }
  }

  // More cleanup.
  int_int_unordered_map_free(node_owners);

  // Exchange and compare the node owners array to make sure everyone gets 
  // the actual owners.
  int* remote_node_owners[num_neighbors];
  for (int p = 0; p < num_neighbors; ++p)
  {
    int num_nodes = matched_node_owners[p]->size;
    remote_node_owners[p] = polymec_malloc(sizeof(int) * num_nodes);
    MPI_Irecv(remote_node_owners[p], num_nodes, MPI_INT, proc, 0, 
              mesh->comm, &requests[p]);
    MPI_Isend(matched_node_owners[p]->data, num_nodes, MPI_INT, proc, 0, 
              mesh->comm, &requests[p + num_neighbors]);
  }
  MPI_Waitall(2 * num_neighbors, requests, statuses);

  for (int p = 0; p < num_neighbors; ++p)
  {
    int num_nodes = matched_node_owners[p]->size;
    for (int n = 0; n < num_nodes; ++n)
    {
      matched_node_owners[p]->data[n] = MIN(matched_node_owners[p]->data[n],
                                            remote_node_owners[p][n]); 
    }
    polymec_free(remote_node_owners[p]);
  }

  //------------------------------------------------------------------------
  // Now matched_nodes[p] has indices of nodes corresponding to those on 
  // neighbor p, and matched_node_owners[p] has arrays containing the actual 
  // owning process for each node in matched_nodes[p]. Next, we make sure 
  // that the owners are sent node requests, and that the exchangers are 
  // constructed accordingly.
  //------------------------------------------------------------------------

  // Gather all the neighbors of our neighbors.
  int_array_t* all_neighbors_of_neighbors = int_array_new();
  {
    int num_neighbors_of_neighbors[num_neighbors];
    for (int p = 0; p < num_neighbors; ++p)
    {
      MPI_Irecv(&num_neighbors_of_neighbors[p], 1, MPI_INT, neighbors[p], 0, 
                mesh->comm, &requests[p]);
      MPI_Isend(&num_neighbors, 1, MPI_INT, neighbors[p], 0, 
          mesh->comm, &requests[p + num_neighbors]);
    }
    MPI_Waitall(2 * num_neighbors, requests, statuses);

    int* neighbors_of_neighbors[num_neighbors];
    for (int p = 0; p < num_neighbors; ++p)
    {
      neighbors_of_neighbors[p] = polymec_malloc(sizeof(int) * num_neighbors_of_neighbors[p]);
      MPI_Irecv(neighbors_of_neighbors[p], num_neighbors_of_neighbors[p], 
          MPI_INT, neighbors[p], 0, mesh->comm, &requests[p]);
      MPI_Isend(neighbors, num_neighbors, MPI_INT, neighbors[p], 0, 
          mesh->comm, &requests[p + num_neighbors]);
    }
    MPI_Waitall(2 * num_neighbors, requests, statuses);

    // Mash them all together.
    for (int p = 0; p < num_neighbors; ++p)
    {
      for (int pp = 0; pp < num_neighbors_of_neighbors[p]; ++pp)
      {
        if (neighbors_of_neighbors[p][pp] != rank)
        {
          int_array_append(all_neighbors_of_neighbors, 
                           neighbors_of_neighbors[p][pp]);
        }
      }
      polymec_free(neighbors_of_neighbors[p]);
    }
  }

  // Now we make a list of all the neighbors of neighbors from whom we 
  // expect data and to whom we will send it.
  {
    int nn_of_n = all_neighbors_of_neighbors->size;
    MPI_Request requests[2 * nn_of_n];
    MPI_Status statuses[2 * nn_of_n];
    int num_nodes_to_send[nn_of_n];
    int_array_t* requested_nodes[nn_of_n];
    for (int p = 0; p < nn_of_n; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(&num_nodes_to_send[p], 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
      requested_nodes[p] = int_array_new();
      for (int j = 0; j < matched_node_owners[p]->size; ++j)
      {
        if (matched_node_owners[p]->data[j] == proc)
          int_array_append(requested_nodes[p], matched_nodes[p]->data[j]);
      }
      MPI_Isend(&requested_nodes[p]->size, 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + nn_of_n]);
    }
    MPI_Waitall(2 * nn_of_n, requests, statuses);

    int_array_t* nodes_to_send[nn_of_n];
    for (int p = 0; p < nn_of_n; ++p)
    {
      nodes_to_send[p] = int_array_new();
      int proc = all_neighbors_of_neighbors->data[p];
      int_array_resize(nodes_to_send[p], num_nodes_to_send[p]);
      MPI_Irecv(nodes_to_send[p]->data, nodes_to_send[p]->size, MPI_INT, 
                proc, 0, mesh->comm, &requests[p]);

      MPI_Isend(requested_nodes[p]->data, requested_nodes[p]->size, 
                MPI_INT, proc, 0, mesh->comm, &requests[p + nn_of_n]);
    }
    MPI_Waitall(2 * nn_of_n, requests, statuses);

    // Set up the exchangers.
    for (int p = 0; p < nn_of_n; ++p)
    {
      if (nodes_to_send[p]->size > 0)
      {
        exchanger_set_send(ex, all_neighbors_of_neighbors->data[p],
                           nodes_to_send[p]->data, nodes_to_send[p]->size,
                           true);
      }
      polymec_free(nodes_to_send[p]);
      if (requested_nodes[p]->size > 0)
      {
        exchanger_set_receive(ex, all_neighbors_of_neighbors->data[p],
                              requested_nodes[p]->data, requested_nodes[p]->size,
                              true);
      }
      int_array_free(requested_nodes[p]);
    }
  }

  // Clean up.
  int_array_free(all_neighbors_of_neighbors);
  for (int p = 0; p < num_neighbors; ++p)
  {
    int_array_free(matched_node_owners[p]);
    int_array_free(matched_nodes[p]);
  }
  exchanger_free(face_ex);
#endif
  return ex;
}

adj_graph_t* graph_from_mesh_cells(mesh_t* mesh)
{
  // Create a graph whose vertices are the mesh's cells.
  int rank, nproc;
  MPI_Comm_size(mesh->comm, &nproc);
  MPI_Comm_rank(mesh->comm, &rank);
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

