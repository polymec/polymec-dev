// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array_utils.h"
#include "core/hilbert.h"
#include "core/kd_tree.h"
#include "core/table.h"
#include "core/unordered_set.h"
#include "geometry/polymesh.h"

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

// Mesh storage stuff.
struct polymesh_storage_t
{
  int cell_face_capacity;
  int face_edge_capacity;
  int face_node_capacity;
  exchanger_t* exchangers[4];
};

// Initializes a new storage mechanism for a polymesh.
static polymesh_storage_t* polymesh_storage_new(MPI_Comm comm)
{
  polymesh_storage_t* storage = polymec_malloc(sizeof(polymesh_storage_t));
  storage->cell_face_capacity = 0;
  storage->face_edge_capacity = 0;
  storage->face_node_capacity = 0;
  memset(storage->exchangers, 0, sizeof(exchanger_t*)*4);
  storage->exchangers[(int)POLYMESH_CELL] = exchanger_new(comm);
  return storage;
}

// Frees the given storage mechanism.
static void polymesh_storage_free(polymesh_storage_t* storage)
{
  for (int cent = 0; cent < 4; ++cent)
  {
    if (storage->exchangers[cent] != NULL)
      release_ref(storage->exchangers[cent]);
  }
  polymec_free(storage);
}

polymesh_t* polymesh_new(MPI_Comm comm, int num_cells, int num_ghost_cells, 
                         int num_faces, int num_nodes)
{
  ASSERT(num_cells >= 0);
  ASSERT(num_ghost_cells >= 0);
  ASSERT(num_faces >= 0);
  ASSERT(num_nodes >= 0);

  polymesh_t* mesh = polymec_malloc(sizeof(polymesh_t));
  mesh->comm = comm;

  // NOTE: We round stored elements up to the nearest power of 2.

  // Allocate cell information.
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;
  mesh->cell_face_offsets = polymec_calloc(num_cells+1, sizeof(int));
  int cell_face_cap = round_to_pow2(12 * num_cells);
  mesh->cell_faces = polymec_calloc(cell_face_cap, sizeof(int));

  // Allocate face information.
  mesh->num_faces = num_faces;

  mesh->face_node_offsets = polymec_calloc(num_faces+1, sizeof(int));
  int face_node_cap = round_to_pow2(6 * num_faces);
  mesh->face_nodes = polymec_calloc(face_node_cap, sizeof(int));

  mesh->face_edge_offsets = polymec_calloc(num_faces+1, sizeof(int));
  mesh->face_edges = NULL;

  mesh->face_cells = polymec_malloc(sizeof(int)*2*num_faces);
  for (int f = 0; f < 2*mesh->num_faces; ++f)
    mesh->face_cells[f] = -1;

  // Allocate edge information.
  mesh->num_edges = 0;
  mesh->edge_nodes = NULL;

  // Allocate node information.
  mesh->num_nodes = num_nodes;
  mesh->nodes = polymec_calloc(num_nodes, sizeof(point_t));

  // Allocate geometric data.
  int total_num_cells = num_cells + num_ghost_cells;
  mesh->cell_volumes = polymec_malloc(sizeof(real_t)*total_num_cells);
  mesh->cell_centers = polymec_malloc(sizeof(point_t)*total_num_cells);
  mesh->face_centers = polymec_malloc(sizeof(point_t)*num_faces);
  mesh->face_areas = polymec_malloc(sizeof(real_t)*num_faces);
  mesh->face_normals = polymec_malloc(sizeof(vector_t)*num_faces);

  // Storage information.
  mesh->storage = polymesh_storage_new(mesh->comm);
  mesh->storage->cell_face_capacity = cell_face_cap;
  mesh->storage->face_node_capacity = face_node_cap;
  mesh->storage->face_edge_capacity = 0;

  // Allocate tagging mechanisms.
  mesh->cell_tags = tagger_new();
  mesh->face_tags = tagger_new();
  mesh->edge_tags = tagger_new();
  mesh->node_tags = tagger_new();

  // Now we create a bogus tag that we can use to store mesh properties.
  int* prop_tag = polymesh_create_tag(mesh->cell_tags, "properties", 1);
  prop_tag[0] = 0;

  return mesh;
}

polymesh_t* polymesh_new_with_cell_type(MPI_Comm comm, int num_cells, 
                                        int num_ghost_cells, int num_faces, 
                                        int num_nodes, int num_faces_per_cell,
                                        int num_nodes_per_face)
{
  polymesh_t* mesh = polymesh_new(comm, num_cells, num_ghost_cells, num_faces, num_nodes);

  // Set up connectivity metadata.
  for (int c = 0; c < mesh->num_cells+1; ++c)
    mesh->cell_face_offsets[c] = num_faces_per_cell*c;
  for (int f = 0; f < mesh->num_faces+1; ++f)
    mesh->face_node_offsets[f] = num_nodes_per_face*f;
  polymesh_reserve_connectivity_storage(mesh);
  return mesh;
}

void polymesh_free(polymesh_t* mesh)
{
  ASSERT(mesh != NULL);

  // Destroy tags.
  tagger_free(mesh->node_tags);
  tagger_free(mesh->edge_tags);
  tagger_free(mesh->face_tags);
  tagger_free(mesh->cell_tags);

  // Destroy storage metadata.
  polymesh_storage_free(mesh->storage);

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

bool polymesh_verify_topology(polymesh_t* mesh, 
                              void (*handler)(const char* format, ...))
{
  // All cells must have at least 4 faces.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (polymesh_cell_num_faces(mesh, c) < 4)
    {
      handler("polymesh_verify_topology: polyhedral cell %d has only %d faces.", 
              c, polymesh_cell_num_faces(mesh, c));
      return false;
    }
  }

  // All faces must have at least 3 nodes/edges.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int ne = polymesh_face_num_edges(mesh, f);
    if (ne == 0)
    {
      handler("polymesh_verify_topology: polygonal face %d has no edges!", f);
      return false;
    }
    if (ne < 3)
    {
      handler("polymesh_verify_topology: polygonal face %d has only %d edges.", f, ne);
      return false;
    }
  }

  // Make sure that all the faces attached to this cell have it in their list.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    int pos = 0, f;
    while (polymesh_cell_next_face(mesh, c, &pos, &f))
    {
      if ((mesh->face_cells[2*f] != c) && (mesh->face_cells[2*f+1] != c))
      {
        handler("polymesh_verify_topology: cell %d has face %d in its list of faces, but that "
                "face does not have that cell in its list of cells.", c, f);
        return false;
      }
    }
  }

  // Now go over all faces and make sure that their cells can see them, too.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int pos = 0, ff;
    bool found_face = false;
    while (polymesh_cell_next_face(mesh, mesh->face_cells[2*f], &pos, &ff))
    {
      if (ff == f) 
      {
        found_face = true;
        break;
      }
    }
    if (!found_face)
    {
      handler("polymesh_verify_topology: face %d has cell %d in its list of cells, but that cell "
              "does not have that face in its list of faces.", f, mesh->face_cells[2*f]);
      return false;
    }

    // Check the cell on the other side, too (but only if its a non-ghost cell).
    if ((mesh->face_cells[2*f+1] != -1) && (mesh->face_cells[2*f+1] < mesh->num_cells))
    {
      pos = 0;
      found_face = false;
      while (polymesh_cell_next_face(mesh, mesh->face_cells[2*f+1], &pos, &ff))
      {
        if (ff == f) 
        {
          found_face = true;
          break;
        }
      }
      if (!found_face)
      {
        handler("polymesh_verify_topology: face %d has cell %d in its list of cells, but that cell "
                "does not have that face in its list of faces.", f, mesh->face_cells[2*f+1]);
        return false;
      }
    }
  }

  return true;
}

polymesh_t* polymesh_clone(polymesh_t* mesh)
{
  polymesh_t* clone = polymesh_new(MPI_COMM_WORLD, mesh->num_cells, mesh->num_ghost_cells,
                           mesh->num_faces, mesh->num_nodes);

  // Connectivity metadata.
  memcpy(clone->cell_face_offsets, mesh->cell_face_offsets, sizeof(int)*(mesh->num_cells+1));
  memcpy(clone->face_node_offsets, mesh->face_node_offsets, sizeof(int)*(mesh->num_faces+1));
  polymesh_reserve_connectivity_storage(clone);

  // Actual connectivity.
  int num_cell_faces = clone->cell_face_offsets[clone->num_cells];
  memcpy(clone->cell_faces, mesh->cell_faces, sizeof(int)*num_cell_faces);
  int num_face_nodes = clone->face_node_offsets[clone->num_faces];
  memcpy(clone->face_nodes, mesh->face_nodes, sizeof(int)*num_face_nodes);
  memcpy(clone->face_cells, mesh->face_cells, sizeof(int)*2*clone->num_faces);

  // Node stuff.
  memcpy(clone->nodes, mesh->nodes, sizeof(point_t)*clone->num_nodes);

  // Construct edges.
  polymesh_construct_edges(clone);

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

int* polymesh_create_tag(tagger_t* tagger, const char* tag, size_t num_indices)
{
  return tagger_create_tag(tagger, tag, num_indices);
}

int* polymesh_tag(tagger_t* tagger, const char* tag, size_t* num_indices)
{
  return tagger_tag(tagger, tag, num_indices);
}

bool polymesh_has_tag(tagger_t* tagger, const char* tag)
{
  return tagger_has_tag(tagger, tag);
}

void polymesh_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag)
{
  tagger_rename_tag(tagger, old_tag, new_tag);
}

void polymesh_delete_tag(tagger_t* tagger, const char* tag)
{
  tagger_delete_tag(tagger, tag);
}

bool polymesh_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, size_t* tag_size)
{
  return tagger_next_tag(tagger, pos, tag_name, tag_indices, tag_size);
}

void polymesh_compute_geometry(polymesh_t* mesh)
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
    while (polymesh_cell_next_oriented_face(mesh, cell, &pos, &face))
    {
      int actual_face = (face >= 0) ? face : ~face;
      ASSERT(actual_face < mesh->num_faces);
      point_t xf = {.x = 0.0, .y = 0.0, .z = 0.0};
      int npos = 0, node;
      while (polymesh_face_next_node(mesh, actual_face, &npos, &node))
      {
        ASSERT(node >= 0);
        ASSERT(node < mesh->num_nodes);
        point_t* xn = &mesh->nodes[node];
        xc.x += xn->x;
        xc.y += xn->y;
        xc.z += xn->z;
        // The first cell on this face computes its center.
        if (cell == mesh->face_cells[2*actual_face])
        {
          xf.x += xn->x;
          xf.y += xn->y;
          xf.z += xn->z;
        }
      }
      if (cell == mesh->face_cells[2*actual_face])
      {
        int nn = polymesh_face_num_nodes(mesh, actual_face);
        xf.x /= nn;
        xf.y /= nn;
        xf.z /= nn;
        mesh->face_centers[actual_face] = xf;
      }
      num_cell_nodes += polymesh_face_num_nodes(mesh, actual_face);
      ++num_cell_faces;
    }
    xc.x /= num_cell_nodes;
    xc.y /= num_cell_nodes;
    xc.z /= num_cell_nodes;
    mesh->cell_centers[cell] = xc;
  }

  // Use the preceding geometry to compute face areas and the 
  // cell's volume.
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    mesh->cell_volumes[cell] = 0.0;
    int pos = 0, face;
    while (polymesh_cell_next_oriented_face(mesh, cell, &pos, &face))
    {
      int actual_face = (face >= 0) ? face : ~face;
      real_t face_area = 0.0;
      vector_t face_normal = {0.0, 0.0, 0.0};
      int epos = 0, edge;
      while (polymesh_face_next_edge(mesh, face, &epos, &edge))
      {
        ASSERT(edge >= 0);
        ASSERT(edge < mesh->num_edges);
        // Construct a tetrahedron whose vertices are the cell center, 
        // the face center, and the two nodes of this edge. The volume 
        // of this tetrahedron contributes to the cell volume.
        vector_t v1, v2, v3, v2xv3;
        point_t* xf = &(mesh->face_centers[actual_face]);
        point_displacement(xf, &(mesh->cell_centers[cell]), &v1);
        point_t xn1 = mesh->nodes[mesh->edge_nodes[2*edge]];
        point_t xn2 = mesh->nodes[mesh->edge_nodes[2*edge+1]];
        point_displacement(xf, &xn1, &v2);
        point_displacement(xf, &xn2, &v3);
        vector_cross(&v2, &v3, &v2xv3);
        real_t tet_volume = ABS(vector_dot(&v1, &v2xv3))/6.0;
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
        point_t* xf = &(mesh->face_centers[actual_face]);
        point_displacement(&(mesh->cell_centers[cell]), xf, &outward);
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

void polymesh_construct_edges(polymesh_t* mesh)
{
  ASSERT(mesh->num_edges == 0);
  ASSERT(mesh->edge_nodes == NULL);

  // Allocate initial face->edge storage.
  if (mesh->storage->face_edge_capacity != mesh->storage->face_node_capacity)
  {
    mesh->storage->face_edge_capacity = mesh->storage->face_node_capacity;
    mesh->face_edges = polymec_realloc(mesh->face_edges, sizeof(int) * mesh->storage->face_edge_capacity);
  }

  // Construct edge information.
  int_table_t* edge_for_nodes = int_table_new();
  memcpy(mesh->face_edge_offsets, mesh->face_node_offsets, sizeof(int) * (mesh->num_faces + 1));
  {
    int num_edges = 0;
    for (int f = 0; f < mesh->num_faces; ++f)
    {
      int offset = mesh->face_edge_offsets[f];
      int num_face_edges = mesh->face_edge_offsets[f+1] - offset;
      for (int e = 0; e < num_face_edges; ++e)
      {
        // Make room.
        if (mesh->storage->face_edge_capacity <= offset+e)
        {
          while (mesh->storage->face_edge_capacity <= offset+e)
            mesh->storage->face_edge_capacity *= 2;
          mesh->face_edges = polymec_realloc(mesh->face_edges, sizeof(int) * mesh->storage->face_edge_capacity);
        }

        int n1 = (int)mesh->face_nodes[offset+e];
        int n2 = (int)mesh->face_nodes[offset+(e+1)%num_face_edges];
        if (!int_table_contains(edge_for_nodes, MIN(n1, n2), MAX(n1, n2)))
        {
          int_table_insert(edge_for_nodes, MIN(n1, n2), MAX(n1, n2), num_edges);
          mesh->face_edges[offset+e] = num_edges;
          ++num_edges;
        }
        else
          mesh->face_edges[offset+e] = *int_table_get(edge_for_nodes, MIN(n1, n2), MAX(n1, n2));
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

void polymesh_reserve_connectivity_storage(polymesh_t* mesh)
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

static size_t polymesh_byte_size(void* obj)
{
  polymesh_t* mesh = obj;
  
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

  // Exchanger-related storage.
  serializer_t* ex_s = exchanger_serializer();
  size_t ex_storage = serializer_size(ex_s, polymesh_exchanger(mesh, POLYMESH_CELL));

  return basic_storage + tag_storage + ex_storage;
}

static void* polymesh_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the number of cells, faces, nodes, and allocate a mesh
  // accordingly.
  int num_cells, num_ghost_cells, num_faces, num_nodes;
  byte_array_read_ints(bytes, 1, &num_cells, offset);
  byte_array_read_ints(bytes, 1, &num_ghost_cells, offset);
  byte_array_read_ints(bytes, 1, &num_faces, offset);
  byte_array_read_ints(bytes, 1, &num_nodes, offset);

  polymesh_t* mesh = polymesh_new(MPI_COMM_WORLD, num_cells, num_ghost_cells,
                                  num_faces, num_nodes);

  // Read all the connectivity metadata.
  byte_array_read_ints(bytes, num_cells+1, mesh->cell_face_offsets, offset);
  byte_array_read_ints(bytes, num_faces+1, mesh->face_node_offsets, offset);

  // Make sure that connectivity storage is sufficient.
  polymesh_reserve_connectivity_storage(mesh);

  // Actual connectivity data.
  int num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
  byte_array_read_ints(bytes, num_cell_faces, mesh->cell_faces, offset);
  int num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
  byte_array_read_ints(bytes, num_face_nodes, mesh->face_nodes, offset);
  byte_array_read_ints(bytes, 2*num_faces, mesh->face_cells, offset);

  // Node stuff.
  byte_array_read_points(bytes, num_nodes, mesh->nodes, offset);

  // Construct edges and compute geometry.
  polymesh_construct_edges(mesh);
  polymesh_compute_geometry(mesh);

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
  polymesh_storage_t* storage = mesh->storage;
  byte_array_read_ints(bytes, 1, &storage->cell_face_capacity, offset);
  byte_array_read_ints(bytes, 1, &storage->face_edge_capacity, offset);
  byte_array_read_ints(bytes, 1, &storage->face_node_capacity, offset);
  release_ref(storage->exchangers[(int)POLYMESH_CELL]);
  ser = exchanger_serializer();
  storage->exchangers[(int)POLYMESH_CELL] = serializer_read(ser, bytes, offset);

  return mesh;
}

static void polymesh_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  polymesh_t* mesh = obj;

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
  polymesh_storage_t* storage = mesh->storage;
  byte_array_write_ints(bytes, 1, &storage->cell_face_capacity, offset);
  byte_array_write_ints(bytes, 1, &storage->face_edge_capacity, offset);
  byte_array_write_ints(bytes, 1, &storage->face_node_capacity, offset);
  ser = exchanger_serializer();
  serializer_write(ser, polymesh_exchanger(mesh, POLYMESH_CELL), bytes, offset);
}

serializer_t* polymesh_serializer()
{
  return serializer_new("mesh", polymesh_byte_size, polymesh_byte_read, polymesh_byte_write, NULL);
}

static exchanger_t* create_face_exchanger(polymesh_t* mesh)
{
  exchanger_t* ex = exchanger_new(mesh->comm);

#if POLYMEC_HAVE_MPI

  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return ex;

  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  // Get the mesh cell exchanger.
  exchanger_t* cell_ex = polymesh_exchanger(mesh, POLYMESH_CELL);

  // Traverse the send/receive cells and create send and receive "faces."
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  int pos = 0, proc, *s_indices, s_size;
  int_unordered_set_t* contributed_to_self = int_unordered_set_new();
  while (exchanger_next_send(cell_ex, &pos, &proc, &s_indices, &s_size))
  {
    // Get the receive transaction corresponding to this send one.
    int *r_indices, r_size; 
    bool have_it = exchanger_get_receive(cell_ex, proc, &r_indices, &r_size);
    if (!have_it)
      polymec_error("create_face_exchanger: Couldn't establish communication!");

    // Now find the face separating each pair of cells.
    ASSERT(r_size == s_size);
    for (int i = 0; i < s_size; ++i)
    {
      int s_cell = s_indices[i];
      int r_cell = r_indices[i];

      int fpos = 0, face;
      while (polymesh_cell_next_face(mesh, s_cell, &fpos, &face))
      {
        int neighbor = polymesh_face_opp_cell(mesh, face, s_cell);
        if (neighbor == r_cell)
        {
          exchanger_proc_map_add_index(send_map, proc, face);
          exchanger_proc_map_add_index(receive_map, proc, face);
          if (!int_unordered_set_contains(contributed_to_self, face))
          {
            exchanger_proc_map_add_index(send_map, rank, face);
            exchanger_proc_map_add_index(receive_map, rank, face);
            int_unordered_set_insert(contributed_to_self, face);
          }
          break;
        }
      }
    }
  }
  exchanger_set_sends(ex, send_map);
  exchanger_set_receives(ex, receive_map);
  int_unordered_set_free(contributed_to_self);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(ex, EXCHANGER_MIN_RANK);
#endif
  return ex;
}

static exchanger_t* create_edge_exchanger(polymesh_t* mesh)
{
  // If we haven't constructed edges yet, do so now.
  if (mesh->num_edges == 0)
    polymesh_construct_edges(mesh);

  exchanger_t* ex = exchanger_new(mesh->comm);
#if POLYMEC_HAVE_MPI
  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return ex;

  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  // The process that owns an edge is the process with minimum MPI rank 
  // associated with a node attached to the edge.
  // Traverse the nodes and associate them with their owning processes.
  // NOTE: We abuse an exchanger_proc_map here to map node indices -> 
  // lists of processes.
#define elem_proc_map_t exchanger_proc_map_t
#define elem_proc_map_new exchanger_proc_map_new
#define elem_proc_map_get exchanger_proc_map_get
#define elem_proc_map_next exchanger_proc_map_next
#define elem_proc_map_add_proc exchanger_proc_map_add_index
#define elem_proc_map_free exchanger_proc_map_free
  exchanger_t* node_ex = polymesh_exchanger(mesh, POLYMESH_NODE);
  elem_proc_map_t* node_procs = elem_proc_map_new();
  {
    int pos = 0, proc, *indices, num_indices;
    while (exchanger_next_receive(node_ex, &pos, &proc, &indices, &num_indices))
    {
      for (int i = 0; i < num_indices; ++i)
        elem_proc_map_add_proc(node_procs, indices[i], proc);
    }
  }

  // Now go over our edges and associate them with the processes on 
  // which the edge lives.
  elem_proc_map_t* edge_procs = elem_proc_map_new();
  int_array_t* edge_indices = int_array_new();
  point_array_t* edge_centers = point_array_new();
  for (int e = 0; e < mesh->num_edges; ++e)
  {
    int n1 = mesh->edge_nodes[2*e];
    int_array_t** n1_procs_p = elem_proc_map_get(node_procs, n1);
    if (n1_procs_p != NULL)
    {
      int_array_t* n1_procs = *n1_procs_p;
      int n2 = mesh->edge_nodes[2*e+1];
      int_array_t** n2_procs_p = elem_proc_map_get(node_procs, n2);
      if (n2_procs_p != NULL)
      {
        int_array_t* n2_procs = *n2_procs_p;

        // The edge must be present on each of the processes we associate
        // with it, so we add the intersection of the two sets of processes
        // n1_procs and n2_procs.
        bool first_entry = true;
        for (size_t i = 0; i < n1_procs->size; ++i)
        {
          int proc = n1_procs->data[i];
          if ((proc != rank) && (int_array_find(n2_procs, proc, int_cmp) != NULL))
          {
            if (first_entry)
            {
              // Jot down the index for this edge and compute its center position.
              int_array_append(edge_indices, e);
              point_t* x1 = &mesh->nodes[n1];
              point_t* x2 = &mesh->nodes[n2];
              point_t xe = {.x = 0.5*(x1->x + x2->x),
                            .y = 0.5*(x1->y + x2->y),
                            .z = 0.5*(x1->z + x2->z)};
              point_array_append(edge_centers, xe);
              first_entry = false;
            }

            elem_proc_map_add_proc(edge_procs, e, proc);
          }
        }
      }
    }
  }
  elem_proc_map_free(node_procs);

  // Sort our edge indices so that they are in Hilbert order.
  {
    bbox_t bbox = {.x1 = REAL_MAX, .x2 = -REAL_MAX,
                   .y1 = REAL_MAX, .y2 = -REAL_MAX,
                   .z1 = REAL_MAX, .z2 = -REAL_MAX};
    for (int i = 0; i < edge_centers->size; ++i)
      bbox_grow(&bbox, &(edge_centers->data[i]));
    hilbert_t* curve = hilbert_new(&bbox);
    hilbert_sort_points(curve, edge_centers->data, edge_indices->data, 
                        edge_centers->size);
    release_ref(curve);
  }

  // Construct the send and receive maps for edges and build the exchanger.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  {
    int_unordered_set_t* contributed_to_self = int_unordered_set_new();
    for (size_t e = 0; e < edge_indices->size; ++e)
    {
      int edge = edge_indices->data[e];
      if (!int_unordered_set_contains(contributed_to_self, edge))
      {
        exchanger_proc_map_add_index(send_map, rank, edge);
        exchanger_proc_map_add_index(receive_map, rank, edge);
        int_unordered_set_insert(contributed_to_self, edge);
      }
      int_array_t* procs = *elem_proc_map_get(edge_procs, edge);
      for (size_t i = 0; i < procs->size; ++i)
      {
        exchanger_proc_map_add_index(send_map, procs->data[i], edge);
        exchanger_proc_map_add_index(receive_map, procs->data[i], edge);
      }
    }
    int_unordered_set_free(contributed_to_self);
  }
  exchanger_set_sends(ex, send_map);
  exchanger_set_receives(ex, receive_map);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(ex, EXCHANGER_MIN_RANK);

  // Clean up.
  point_array_free(edge_centers);
  int_array_free(edge_indices);
  elem_proc_map_free(edge_procs);
#undef node_proc_map_t
#undef node_proc_map_new 
#undef node_proc_map_get 
#undef node_proc_map_add_proc 
#undef node_proc_map_free 

#endif
  return ex;
}

static exchanger_t* create_node_exchanger(polymesh_t* mesh)
{
  exchanger_t* ex = exchanger_new(mesh->comm);
#if POLYMEC_HAVE_MPI

  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return ex;

  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  // Generate a list of all the neighbors of our set of neighbors.
  int_array_t* all_neighbors_of_neighbors = int_array_new();
  exchanger_t* face_ex = polymesh_exchanger(mesh, POLYMESH_FACE);
  int num_neighbors = exchanger_num_sends(face_ex);
  int neighbors[num_neighbors];
  {
    MPI_Request requests[2*num_neighbors];
    MPI_Status statuses[2*num_neighbors];

    // Identify our own neighbors.
    int pos = 0, proc, *indices, size, p = 0;
    while (exchanger_next_send(face_ex, &pos, &proc, &indices, &size))
      neighbors[p++] = proc;

    // Get the number of neighbors for our pth neighbor.
    int num_neighbors_of_neighbors[num_neighbors];
    for (int pp = 0; pp < num_neighbors; ++pp)
    {
      MPI_Irecv(&num_neighbors_of_neighbors[pp], 1, MPI_INT, neighbors[pp], 0, 
                mesh->comm, &requests[pp]);
      MPI_Isend(&num_neighbors, 1, MPI_INT, neighbors[pp], 0, 
          mesh->comm, &requests[pp + num_neighbors]);
    }
    MPI_Waitall(2 * num_neighbors, requests, statuses);

    // Get the ranks of the neighbors for our pth neighbor.
    int* neighbors_of_neighbors[num_neighbors];
    for (int pp = 0; pp < num_neighbors; ++pp)
    {
      neighbors_of_neighbors[pp] = polymec_malloc(sizeof(int) * num_neighbors_of_neighbors[pp]);
      MPI_Irecv(neighbors_of_neighbors[pp], num_neighbors_of_neighbors[pp], 
          MPI_INT, neighbors[pp], 0, mesh->comm, &requests[pp]);
      MPI_Isend(neighbors, num_neighbors, MPI_INT, neighbors[pp], 0, 
          mesh->comm, &requests[pp + num_neighbors]);
    }
    MPI_Waitall(2 * num_neighbors, requests, statuses);

    // Mash them all together.
    int_unordered_set_t* neighbor_neighbor_set = int_unordered_set_new();
    for (int pp = 0; pp < num_neighbors; ++pp)
    {
      int_unordered_set_insert(neighbor_neighbor_set, neighbors[pp]);
      for (int ppp = 0; ppp < num_neighbors_of_neighbors[pp]; ++ppp)
      {
        if (neighbors_of_neighbors[pp][ppp] != rank)
        {
          int_unordered_set_insert(neighbor_neighbor_set, 
                                   neighbors_of_neighbors[pp][ppp]);
        }
      }
      polymec_free(neighbors_of_neighbors[pp]);
    }
    int npos = 0, neighbor;
    while (int_unordered_set_next(neighbor_neighbor_set, &npos, &neighbor))
      int_array_append(all_neighbors_of_neighbors, neighbor);
    int_unordered_set_free(neighbor_neighbor_set);
  }

  // Make a list of all the nodes that can be communicated to neighboring 
  // processes.
  int_array_t* my_node_indices = int_array_new();
  point_array_t* my_nodes = point_array_new();
  {
    int_unordered_set_t* node_set = int_unordered_set_new();
    int pos = 0, proc, *indices, size;
    while (exchanger_next_send(face_ex, &pos, &proc, &indices, &size))
    {
      for (int f = 0; f < size; ++f)
      {
        int face = indices[f];
        int npos = 0, node;
        while (polymesh_face_next_node(mesh, face, &npos, &node))
          int_unordered_set_insert(node_set, node);
      }
    }
    int npos = 0, node;
    while (int_unordered_set_next(node_set, &npos, &node))
    {
      int_array_append(my_node_indices, node);
      point_array_append(my_nodes, mesh->nodes[node]);
    }
    int_unordered_set_free(node_set);
  }

  // Sort our nodes so that they are in Hilbert order.
  {
    bbox_t bbox = {.x1 = REAL_MAX, .x2 = -REAL_MAX,
                   .y1 = REAL_MAX, .y2 = -REAL_MAX,
                   .z1 = REAL_MAX, .z2 = -REAL_MAX};
    for (int i = 0; i < my_nodes->size; ++i)
      bbox_grow(&bbox, &(my_nodes->data[i]));
    hilbert_t* curve = hilbert_new(&bbox);
    hilbert_sort_points(curve, my_nodes->data, my_node_indices->data, 
                        my_nodes->size);
    release_ref(curve);
  }

  // Now send/receive the positions of all nodes that can interact with 
  // neighbors of our neighbors. 
  int num_neighbor_neighbors = (int)all_neighbors_of_neighbors->size;
  MPI_Request requests[2*num_neighbor_neighbors];
  MPI_Status statuses[2*num_neighbor_neighbors];
  int_ptr_unordered_map_t* neighbor_neighbor_nodes = int_ptr_unordered_map_new();
  {
    // How many nodes does each neighbor neighbor have?
    int num_neighbor_neighbor_nodes[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(&num_neighbor_neighbor_nodes[p], 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
    }
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Isend(&(my_nodes->size), 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);

    // Get 'em!
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      point_array_t* nn_nodes = point_array_new();
      point_array_resize(nn_nodes, num_neighbor_neighbor_nodes[p]);
      int_ptr_unordered_map_insert_with_v_dtor(neighbor_neighbor_nodes, proc, nn_nodes, DTOR(point_array_free));
      MPI_Irecv(nn_nodes->data, 3*num_neighbor_neighbor_nodes[p], MPI_REAL_T, proc, 0, 
                mesh->comm, &requests[p]);
    }

    // Send 'em!
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Isend(my_nodes->data, (int)(3*my_nodes->size), MPI_REAL_T, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);
  }

  // At this point, neighbor_neighbor_nodes maps the ranks of all processes
  // we can possibly interact with to the positions of the nodes on those 
  // processes. 

  // Now we cut out the nodes that don't match up.
  {
    // Post the receives for the numbers of culled nodes.
    int num_culled_nodes[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(&(num_culled_nodes[p]), 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
    }

    // Figure out which nodes we will cull and send the number.
    // culled_nodes[p] contains a list of the nodes on neighbor p 
    // that DON'T correspond to any of my_nodes.
    int_array_t* culled_nodes[num_neighbor_neighbors];
    real_t tolerance = 1e-8; // FIXME: This is bad.
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      point_array_t* nn_nodes = *int_ptr_unordered_map_get(neighbor_neighbor_nodes, proc);
      culled_nodes[p] = int_array_new();
      int_unordered_set_t* kept_nodes = int_unordered_set_new();
      for (int i = 0; i < my_nodes->size; ++i)
      {
        point_t* xi = &my_nodes->data[i];
        for (int j = 0; j < nn_nodes->size; ++j)
        {
          point_t* xj = &nn_nodes->data[j];
          if (point_distance(xi, xj) <= tolerance)
            int_unordered_set_insert(kept_nodes, j);
        }
      }
      for (int i = 0; i < nn_nodes->size; ++i)
      {
        if (!int_unordered_set_contains(kept_nodes, i))
          int_array_append(culled_nodes[p], i);
      }
      int_unordered_set_free(kept_nodes);
      MPI_Isend(&culled_nodes[p]->size, 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);

    // Now receive/send the culled nodes.
    // my_culled_nodes[p] is the list of indices of nodes in my_nodes
    // that DON'T correspond to any of the nodes on neighbor p.
    int_array_t* my_culled_nodes[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      my_culled_nodes[p] = int_array_new();
      int_array_resize(my_culled_nodes[p], num_culled_nodes[p]);
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(my_culled_nodes[p]->data, num_culled_nodes[p], MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
    }
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Isend(culled_nodes[p]->data, (int)culled_nodes[p]->size, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);

    // Organized the culled nodes into sets for querying.
    int_unordered_set_t* my_culled_node_sets[num_neighbor_neighbors];
    int_unordered_set_t* their_culled_node_sets[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      my_culled_node_sets[p] = int_unordered_set_new();
      for (int i = 0; i < my_culled_nodes[p]->size; ++i)
        int_unordered_set_insert(my_culled_node_sets[p], my_culled_nodes[p]->data[i]);
      int_array_free(my_culled_nodes[p]);
      their_culled_node_sets[p] = int_unordered_set_new();
      for (int i = 0; i < culled_nodes[p]->size; ++i)
        int_unordered_set_insert(their_culled_node_sets[p], culled_nodes[p]->data[i]);
      int_array_free(culled_nodes[p]);
    }

    // Create a kd tree that stores all the nodes in 3D space.
    kd_tree_t* node_tree = kd_tree_new(mesh->nodes, mesh->num_nodes);

    // Now set up the exchanger send/receive maps.
    exchanger_proc_map_t* send_map = exchanger_proc_map_new();
    exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
    int_unordered_set_t* contributed_to_self = int_unordered_set_new();
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];

      // Set up send nodes:
      for (int i = 0; i < my_nodes->size; ++i)
      {
        int node = my_node_indices->data[i];

        // Set up the self contribution.
        if (!int_unordered_set_contains(contributed_to_self, node))
        {
          exchanger_proc_map_add_index(send_map, rank, node);
          exchanger_proc_map_add_index(receive_map, rank, node);
          int_unordered_set_insert(contributed_to_self, node);
        }

        // Get other contributions.
        if (!int_unordered_set_contains(my_culled_node_sets[p], i))
          exchanger_proc_map_add_index(send_map, proc, node);
      }
      int_unordered_set_free(my_culled_node_sets[p]);

      // Set up receive nodes.
      point_array_t* their_nodes = *int_ptr_unordered_map_get(neighbor_neighbor_nodes, proc);
      int sorted_indices[their_nodes->size];
      for (int i = 0; i < their_nodes->size; ++i)
        sorted_indices[i] = i;
      int_qsort(sorted_indices, their_nodes->size);

      for (int i = 0; i < their_nodes->size; ++i)
      {
        int j = sorted_indices[i];
        if (!int_unordered_set_contains(their_culled_node_sets[p], j))
        {
          // Find the node in "their_nodes" that matches our local node and 
          // add it to our list of receive nodes.
          int node = kd_tree_nearest(node_tree, &their_nodes->data[j]);
          exchanger_proc_map_add_index(receive_map, proc, node);
        }
      }

      int_unordered_set_free(their_culled_node_sets[p]);
    }
    int_unordered_set_free(contributed_to_self);

    kd_tree_free(node_tree);
    int_array_free(my_node_indices);
    point_array_free(my_nodes);

    // Add the maps to the exchanger.
    exchanger_set_sends(ex, send_map);
    exchanger_set_receives(ex, receive_map);
  }

  // Clean up.
  int_ptr_unordered_map_free(neighbor_neighbor_nodes);
  int_array_free(all_neighbors_of_neighbors);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(ex, EXCHANGER_MIN_RANK);
#endif
  return ex;
}

adj_graph_t* graph_from_polymesh_cells(polymesh_t* mesh)
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

exchanger_t* polymesh_exchanger(polymesh_t* mesh,
                                polymesh_centering_t centering)
{
  if ((centering == POLYMESH_FACE) &&
      (mesh->storage->exchangers[(int)centering] == NULL))
    mesh->storage->exchangers[(int)centering] = create_face_exchanger(mesh);
  else if ((centering == POLYMESH_EDGE) &&
           (mesh->storage->exchangers[(int)centering] == NULL))
    mesh->storage->exchangers[(int)centering] = create_edge_exchanger(mesh);
  else if ((centering == POLYMESH_NODE) && 
           (mesh->storage->exchangers[(int)centering] == NULL))
    mesh->storage->exchangers[(int)centering] = create_node_exchanger(mesh);
  return mesh->storage->exchangers[(int)centering];
}

void polymesh_set_exchanger(polymesh_t* mesh, 
                            polymesh_centering_t centering, 
                            exchanger_t* ex);
void polymesh_set_exchanger(polymesh_t* mesh, 
                            polymesh_centering_t centering, 
                            exchanger_t* ex)
{
  ASSERT(ex != NULL);
  release_ref(mesh->storage->exchangers[(int)centering]);
  mesh->storage->exchangers[(int)centering] = ex;
}

