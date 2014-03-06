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

// This implementation of the Voronoi tessellator uses Tetgen. It is not built 
// unless a Tetgen tarball was found in the 3rdparty/ directory.
// Please see the license file in the Tetgen tarball for license information 
// on Tetgen.

#include "core/mesh_storage.h"
#include "core/table.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/text_file_buffer.h"
#include "geometry/create_tetgen_mesh.h"

typedef struct
{
  int num_nodes; // 4 for order 1, 10 for order 2.
  int nodes[10];  
  int attribute; // 0 for none, positive for actual attribute.
  int neighbors[4]; // Neighboring tetrahedra.
} tet_t;

typedef struct
{
  int num_nodes; // 3 for order 1, 6 for order 2.
  int nodes[6];  
  int boundary_marker; // 0 for none, positive for actual attribute.
} tet_face_t;

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

static point_t* read_nodes(const char* node_file, int* num_nodes)
{
  *num_nodes = -1;
  point_t* nodes = NULL;

  text_file_buffer_t* buffer = text_file_buffer_new(node_file);
  if (buffer == NULL)
    polymec_error("TetGen node file '%s' not found.", node_file);
  int nodes_read = 0, pos = 0, line_length;
  char* line;
  while (text_file_buffer_next(buffer, &pos, &line, &line_length))
  {
    // Skip lines starting with #.
    if (line[0] == '#') continue;

    // Look for the header if we haven't already read it.
    if (*num_nodes == -1)
    {
      int dim, num_attributes, num_boundary_markers; // ignored
      int num_items = sscanf(line, "%d %d %d %d\n", num_nodes, &dim, 
                             &num_attributes, &num_boundary_markers);
      if (num_items != 4)
        polymec_error("Node file has bad header.");
      if (dim != 3)
        polymec_error("Node file is not 3-dimensional.");
      nodes = malloc(sizeof(point_t) * (*num_nodes));
      continue;
    }

    // Read the coordinates of the next node
    char garbage[4096];
    int node_id;
    int num_items = sscanf(line, "%d %lg %lg %lg%s\n", &node_id, 
                           &nodes[nodes_read].x, &nodes[nodes_read].y, 
                           &nodes[nodes_read].z, garbage);
    if (node_id != (nodes_read+1))
      polymec_error("Bad node ID after %d nodes read: %d.\n", nodes_read, node_id);
    if (num_items < 4)
      polymec_error("Bad line in nodes file after %d nodes read.\n", nodes_read);
    if (*num_nodes <= 0)
      polymec_error("Bad number of nodes in node file: %d.", *num_nodes);

    ++nodes_read;
    if (nodes_read == *num_nodes) break;
  }
  text_file_buffer_free(buffer);

  if (nodes_read != *num_nodes)
    polymec_error("Node file claims to contain %d nodes, but %d were read.", *num_nodes, nodes_read);
  return nodes;
}

static tet_t* read_tets(const char* tet_file, int* num_tets)
{
  *num_tets = -1;
  tet_t* tets = NULL;
  text_file_buffer_t* buffer = text_file_buffer_new(tet_file);
  if (buffer == NULL)
    polymec_error("TetGen element file '%s' not found.", tet_file);
  int tets_read = 0, nodes_per_tet = 4, region_attribute = 0, pos = 0, line_length;
  char* line;
  while (text_file_buffer_next(buffer, &pos, &line, &line_length))
  {
    // Skip lines starting with #.
    if (line[0] == '#') continue;

    // Look for the header if we haven't already read it.
    if (*num_tets == -1)
    {
      int num_items = sscanf(line, "%d %d %d\n", num_tets, 
                             &nodes_per_tet, &region_attribute);
      if (num_items != 3)
        polymec_error("Element file has bad header.");
      if (*num_tets <= 0)
        polymec_error("Bad number of tets in element file: %d.", *num_tets);
      if ((nodes_per_tet != 4) && (nodes_per_tet != 10))
        polymec_error("Bad number of nodes per tet: %d (must be 4 or 10).", nodes_per_tet);
      tets = malloc(sizeof(tet_t) * (*num_tets));
      continue;
    }

    // Read the node indices of the next tet.
    tet_t* tet = &tets[tets_read];
    char attr_str[128];
    int tet_id;
    if (nodes_per_tet == 4)
    {
      tet->num_nodes = 4;
      int num_items = sscanf(line, "%d %d %d %d %d%s\n", &tet_id, 
                             &tet->nodes[0], &tet->nodes[1], 
                             &tet->nodes[2], &tet->nodes[3],
                             attr_str);
      if (num_items < 5)
        polymec_error("Bad line in element file after %d tets read.\n", tets_read);
      if (num_items == 6)
        tet->attribute = atoi(attr_str);
    }
    else
    {
      tet->num_nodes = 10;
      int num_items = sscanf(line, "%d %d %d %d %d %d %d %d %d %d %d%s\n", 
                             &tet_id, &tet->nodes[0], &tet->nodes[1], 
                             &tet->nodes[2], &tet->nodes[3], &tet->nodes[4],
                             &tet->nodes[5], &tet->nodes[6], &tet->nodes[7],
                             &tet->nodes[8], &tet->nodes[9], attr_str);
      if (num_items < 11)
        polymec_error("Bad line in element file after %d tets read.\n", tets_read);
      if (num_items == 12)
        tet->attribute = atoi(attr_str);
    }
    if (tet_id != (tets_read+1))
      polymec_error("Bad tet ID after %d tets read: %d.\n", tets_read, tet_id);
    ++tets_read;
    if (tets_read == *num_tets) break;
  }
  text_file_buffer_free(buffer);
  if (tets_read != *num_tets)
    polymec_error("Element file claims to contain %d tets, but %d were read.", *num_tets, tets_read);

  // TetGen's indices are 1-based, so correct them.
  for (int t = 0; t < *num_tets; ++t)
  {
    tet_t* tet = &tets[t];
    for (int n = 0; n < tet->num_nodes; ++n)
      tet->nodes[n] -= 1;
  }

  return tets;
}

static tet_face_t* read_faces(const char* face_file, int nodes_per_face, int* num_faces)
{
  *num_faces = -1;
  tet_face_t* faces = NULL;
  text_file_buffer_t* buffer = text_file_buffer_new(face_file);
  if (buffer == NULL)
    polymec_error("TetGen face file '%s' not found.", face_file);
  int faces_read = 0, boundary_marker = 0, pos = 0, line_length;
  char* line;
  while (text_file_buffer_next(buffer, &pos, &line, &line_length))
  {
    // Skip lines starting with #.
    if (line[0] == '#') continue;

    // Look for the header if we haven't already read it.
    if (*num_faces == -1)
    {
      int num_items = sscanf(line, "%d %d\n", num_faces, &boundary_marker);
      if (num_items != 2)
        polymec_error("Face file has bad header.");
      if (*num_faces <= 0)
        polymec_error("Bad number of faces in face file: %d.", *num_faces);
      faces = malloc(sizeof(tet_face_t) * (*num_faces));
      continue;
    }

    // Read the node indices of the next face.
    tet_face_t* face = &faces[faces_read];
    char marker_str[128];
    int face_id;
    if (nodes_per_face == 3)
    {
      face->num_nodes = 3;
      int num_items = sscanf(line, "%d %d %d %d%s\n", &face_id, 
                             &face->nodes[0], &face->nodes[1], 
                             &face->nodes[2], marker_str);
      if (num_items < 4)
        polymec_error("Bad line in face file after %d faces read.\n", faces_read);
      if (num_items == 5)
        face->boundary_marker = atoi(marker_str);
    }
    else
    {
      face->num_nodes = 6;
      int num_items = sscanf(line, "%d %d %d %d %d %d %d%s\n", 
                             &face_id, &face->nodes[0], &face->nodes[1], 
                             &face->nodes[2], &face->nodes[3], &face->nodes[4],
                             &face->nodes[5], marker_str);
      if (num_items < 7)
        polymec_error("Bad line in face file after %d faces read.\n", faces_read);
      if (num_items == 8)
        face->boundary_marker = atoi(marker_str);
    }
    if (face_id != (faces_read+1))
      polymec_error("Bad face ID after %d faces read: %d.\n", faces_read, face_id);
    ++faces_read;
    if (faces_read == *num_faces) break;
  }
  text_file_buffer_free(buffer);
  if (faces_read != *num_faces)
    polymec_error("Face file claims to contain %d faces, but %d were read.", *num_faces, faces_read);

  // TetGen's indices are 1-based, so correct them.
  for (int f = 0; f < *num_faces; ++f)
  {
    tet_face_t* face = &faces[f];
    for (int n = 0; n < face->num_nodes; ++n)
      face->nodes[n] -= 1;
  }

  return faces;
}

static void read_neighbors(const char* neigh_file, tet_t* tets, int num_tets)
{
  text_file_buffer_t* buffer = text_file_buffer_new(neigh_file);
  if (buffer == NULL)
    polymec_error("TetGen neighbor file '%s' not found.", neigh_file);
  int tets_read = 0, pos = 0, line_length, num_entries = -1;
  char* line;
  while (text_file_buffer_next(buffer, &pos, &line, &line_length))
  {
    // Skip lines starting with #.
    if (line[0] == '#') continue;

    // Look for the header if we haven't already read it.
    if (num_entries == -1)
    {
      int four;
      int num_items = sscanf(line, "%d %d\n", &num_entries, &four);
      if (num_items != 2)
        polymec_error("Face file has bad header.");
      if (num_entries != num_tets)
        polymec_error("Number of neighbor entries (%d) in neigh file does not match number of tets (%d).", num_entries, num_tets);
      if (four != 4)
        polymec_error("Second value in header must be 4.");
      continue;
    }

    // Read the neighbor indices of the next tet.
    tet_t* tet = &tets[tets_read];
    int tet_id;
    char junk_str[1024];
    int num_items = sscanf(line, "%d %d %d %d %d%s\n", &tet_id, 
                           &tet->neighbors[0], &tet->neighbors[1], 
                           &tet->neighbors[2], &tet->neighbors[3], 
                           junk_str);
    if (num_items < 5)
        polymec_error("Bad line in neighbors file after %d tets read.\n", tets_read);
    if (tet_id != (tets_read+1))
      polymec_error("Bad tet ID after %d tet read: %d.\n", tets_read, tet_id);
    ++tets_read;
    if (tets_read == num_tets) break;
  }
  text_file_buffer_free(buffer);
  if (tets_read != num_tets)
    polymec_error("Neighbor file has %d tets, but needs %d.", tets_read, num_tets);

  // TetGen's indices are 1-based, so correct them.
  for (int t = 0; t < num_tets; ++t)
  {
    tet_t* tet = &tets[t];
    for (int n = 0; n < 4; ++n)
      tet->neighbors[n] -= 1;
  }
}

static int_table_t* gather_edges(tet_face_t* faces, 
                                 int num_faces,
                                 int* num_edges)
{
  int_table_t* edge_for_nodes = int_table_new();
  *num_edges = 0;
  for (int f = 0; f < num_faces; ++f)
  {
    tet_face_t* face = &faces[f];
    for (int e = 0; e < face->num_nodes; ++e)
    {
      int n1 = face->nodes[e];
      int n2 = face->nodes[(e+1)%face->num_nodes];
      int_table_insert(edge_for_nodes, n1, n2, *num_edges);
      *num_edges += 1;
    }
  }
  return edge_for_nodes;
}

static bool even_permutation(int n1, int n2, int n3, int m1, int m2, int m3)
{
  return (((n1 == m1) && (n2 == m2) && (n3 == m3)) || 
          ((n1 == m2) && (n2 == m3) && (n3 == m1)) ||
          ((n1 == m3) && (n2 == m1) && (n3 == m2)));
}

static bool odd_permutation(int n1, int n2, int n3, int m1, int m2, int m3)
{
  return (((n1 == m2) && (n2 == m1) && (n3 == m3)) || 
          ((n1 == m1) && (n2 == m3) && (n3 == m2)) ||
          ((n1 == m3) && (n2 == m2) && (n3 == m1)));
}

static bool nodes_are_even_permutation_of_tet_face(int n1, int n2, int n3, tet_t* tet)
{
  // According to TetGen's indexing scheme, a tet has 4 faces with 
  // the following locally-indexed nodes:
  // 1. 0, 1, 3
  // 2. 1, 2, 3
  // 3. 0, 2, 1
  // 4. 2, 0, 3
  return (even_permutation(n1, n2, n3, tet->nodes[0], tet->nodes[1], tet->nodes[3]) ||
          even_permutation(n1, n2, n3, tet->nodes[1], tet->nodes[2], tet->nodes[3]) ||
          even_permutation(n1, n2, n3, tet->nodes[0], tet->nodes[2], tet->nodes[1]) ||
          even_permutation(n1, n2, n3, tet->nodes[2], tet->nodes[0], tet->nodes[3]));
}

static int find_face_with_nodes(tet_face_t* faces, int num_faces, int n1, int n2, int n3)
{
  for (int f = 0; f < num_faces; ++f)
  {
    tet_face_t* face = &faces[f];
    if (even_permutation(face->nodes[0], face->nodes[1], face->nodes[2], n1, n2, n3) ||
        odd_permutation(face->nodes[0], face->nodes[1], face->nodes[2], n1, n2, n3))
      return f;
  }
  return -1;
}

mesh_t* create_tetgen_mesh(MPI_Comm comm, 
                           const char* node_file,
                           const char* ele_file,
                           const char* face_file,
                           const char* neigh_file)
{
  int num_nodes;
  point_t* nodes = read_nodes(node_file, &num_nodes);

  int num_tets;
  tet_t* tets = read_tets(ele_file, &num_tets);

  int num_faces;
  int nodes_per_face = (tets[0].num_nodes == 4) ? 3 : 6;
  tet_face_t* faces = read_faces(face_file, nodes_per_face, &num_faces);

  read_neighbors(neigh_file, tets, num_tets);

  // Compute the number of edges, which we aren't given.
  int num_edges = 0;
  int_table_t* edge_for_nodes = gather_edges(faces, num_faces, &num_edges);

  // Create the mesh.
  mesh_t* mesh = mesh_new(comm, num_tets, 0, num_faces, num_edges, num_nodes);
  
  // Copy node coordinates.
  memcpy(mesh->nodes, nodes, sizeof(point_t) * num_nodes);

  // Edge <-> node connectivity.
  {
    int_table_cell_pos_t pos = int_table_start(edge_for_nodes);
    int n1, n2, e;
    while (int_table_next_cell(edge_for_nodes, &pos, &n1, &n2, &e))
    {
      mesh->edge_nodes[2*e] = n1;
      mesh->edge_nodes[2*e+1] = n2;
    }
  }

  // Face <-> node connectivity.
  mesh->face_node_offsets[0] = 0;
  for (int f = 0; f < num_faces; ++f)
    mesh->face_node_offsets[f+1] = (f+1)*nodes_per_face;
  mesh->face_nodes = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int) * num_faces * nodes_per_face, 0);
  for (int f = 0; f < num_faces; ++f)
  {
    for (int n = 0; n < nodes_per_face; ++n)
      mesh->face_nodes[nodes_per_face*f+n] = faces->nodes[n];
  }
  mesh->storage->face_node_capacity = num_faces * nodes_per_face;

  // Face <-> edge connectivity.
  memset(mesh->face_edge_offsets, 0, sizeof(int) * (mesh->num_faces + 1));
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int Ne = mesh->face_node_offsets[f+1] - mesh->face_node_offsets[f];
    mesh->face_edge_offsets[f+1] = Ne;
    for (int e = 0; e < Ne; ++e)
    {
      int offset = mesh->face_node_offsets[f];
      int n1 = (int)mesh->face_nodes[offset+e];
      int n2 = (int)mesh->face_nodes[offset+(e+1)%Ne];
      int edge_id = *int_table_get(edge_for_nodes, n1, n2);
      mesh->storage->face_edge_capacity = round_to_pow2(edge_id+1);
      mesh->face_edges = ARENA_REALLOC(mesh->arena, mesh->face_edges, sizeof(int) * mesh->storage->face_edge_capacity, 0);
      mesh->face_edges[mesh->face_node_offsets[f+1]+e] = edge_id;
    }
  }

  // Cell <-> face connectivity.
  mesh->cell_face_offsets[0] = 0;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    mesh->cell_face_offsets[c+1] = 4*(c+1);
    for (int f = mesh->cell_face_offsets[c]; f < mesh->cell_face_offsets[c+1]; ++f)
      mesh->cell_faces[f] = -1;
  }
  for (int f = 0; f < num_faces; ++f)
  {
    mesh->face_cells[2*f]   = -1;
    mesh->face_cells[2*f+1] = -1;
  }
  mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int) * 4 * mesh->num_cells, 0);
  {
    // Loop over cells and find the faces connecting them to their neighbors.
    for (int c = 0; c < mesh->num_cells; ++c)
    {
      tet_t* t = &tets[c];

      // Tet face-node table.
      static int tet_face_nodes[4][3] = {{1, 2, 3},  // face 1 has nodes 2, 3, 4
                                         {2, 0, 3},  // face 2 has nodes 3, 1, 4
                                         {0, 1, 3},  // face 3 has nodes 1, 2, 4
                                         {0, 1, 2}}; // face 4 has nodes 1, 2, 3

      // Keep track of faces we've processed.
      int_unordered_set_t* faces_processed = int_unordered_set_new();

      // Figure out each of the connections by examining their common nodes.
      // We use TetGen's indexing scheme (see TetGen documentation), which 
      // states that neighbor n of a tet shares the face of that tet that 
      // is opposite of node n in the tet.
      for (int n = 0; n < 4; ++n)
      {
        // Nodes of cell c on this face.
        int c_n1 = t->nodes[tet_face_nodes[n][0]];
        int c_n2 = t->nodes[tet_face_nodes[n][1]];
        int c_n3 = t->nodes[tet_face_nodes[n][2]];

        // Get the neighbor tet.
        int cn = t->neighbors[n];
        if (cn == -1)
        {
          // This is a boundary face! Find the face and set it up.
          int face = find_face_with_nodes(faces, num_faces, c_n1, c_n2, c_n3);
          mesh->cell_faces[n] = face; // FIXME: or ~face?
          mesh->face_cells[2*face] = c;
          int_unordered_set_insert(faces_processed, face);
        }
        else if (cn > c)
        {
          tet_t* tn = &tets[cn];

          // Find the neighbor index of c within cn.
          int n1 = (tn->neighbors[0] == c) ? 0 :
            (tn->neighbors[1] == c) ? 1 : 
            (tn->neighbors[2] == c) ? 2 : 3;

          // The nodes in these faces should match up.
          ASSERT(t0->nodes[tet_face_nodes[n][0]] == tn->nodes[tet_face_nodes[n1][0]]);
          ASSERT(t0->nodes[tet_face_nodes[n][1]] == tn->nodes[tet_face_nodes[n1][1]]);
          ASSERT(t0->nodes[tet_face_nodes[n][2]] == tn->nodes[tet_face_nodes[n1][2]]);
               
          // Find the face that possesses these 3 nodes.
          int face = find_face_with_nodes(faces, num_faces, c_n1, c_n2, c_n3);

          // Associate the face with both of these cells.
          mesh->cell_faces[n] = face; // FIXME: or ~face?
          mesh->cell_faces[n1] = ~face; // FIXME: or face?

          // Associate the cells with the face.
          mesh->face_cells[2*face]   = c;
          mesh->face_cells[2*face+1] = cn;

          int_unordered_set_insert(faces_processed, face);
        }
      }

      // Clean up.
      int_unordered_set_free(faces_processed);
    }

#if 0
    // Construct a mapping from nodes to sets of cells.
    int_ptr_unordered_map_t* cells_for_node = int_ptr_unordered_map_new();
    for (int c = 0; c < num_tets; ++c)
    {
      tet_t* tet = &tets[c];
      for (int n = 0; n < tet->num_nodes; ++n)
      {
        int node = tet->nodes[n];
        int_unordered_set_t** cells = (int_unordered_set_t**)int_ptr_unordered_map_get(cells_for_node, node);
        if (cells == NULL)
        {
          int_unordered_set_t* cells_for_this_node = int_unordered_set_new();
          int_unordered_set_insert(cells_for_this_node, c);
          int_ptr_unordered_map_insert_with_v_dtor(cells_for_node, node, cells_for_this_node, DTOR(int_unordered_set_free));
        }
        else
        {
          int_unordered_set_t* cells_for_this_node = *cells;
          int_unordered_set_insert(cells_for_this_node, c);
        }
      }
    }

    // Now go over faces and find those cells that possess all of the 
    // nodes in each face. Those cells are the cells that share the face.
    int_unordered_set_t* intersect01 = int_unordered_set_new();
    int_unordered_set_t* intersect012 = int_unordered_set_new();
    for (int f = 0; f < num_faces; ++f)
    {
      tet_face_t* face = &faces[f];
      int n0 = face->nodes[0];
      int_unordered_set_t** node_cells0 = (int_unordered_set_t**)int_ptr_unordered_map_get(cells_for_node, n0);
      ASSERT(node_cells0 != NULL);
      int_unordered_set_t* cells_for_0 = *node_cells0;
      int n1 = face->nodes[1];
      int_unordered_set_t** node_cells1 = (int_unordered_set_t**)int_ptr_unordered_map_get(cells_for_node, n1);
      ASSERT(node_cells1 != NULL);
      int_unordered_set_t* cells_for_1 = *node_cells1;
      int n2 = face->nodes[2];
      int_unordered_set_t** node_cells2 = (int_unordered_set_t**)int_ptr_unordered_map_get(cells_for_node, n2);
      ASSERT(node_cells2 != NULL);
      int_unordered_set_t* cells_for_2 = *node_cells2;
      int_unordered_set_intersection(cells_for_0, cells_for_1, intersect01);
      int_unordered_set_intersection(intersect01, cells_for_2, intersect012);
      ASSERT((intersect012->size == 1) || (intersect012->size == 2));
      int pos = 0, cell, i = 0;
      while (int_unordered_set_next(intersect012, &pos, &cell))
      {
        tet_t* tet = &tets[cell];

        // Hook the cell up to the face.
        mesh->face_cells[2*f+i] = cell;
        ++i;

        // Hook the face up to the cell.
        int j = mesh->cell_face_offsets[cell];
        while (j < mesh->cell_face_offsets[cell+1])
        {
          if (mesh->cell_faces[j] == -1)
          {
            // We have to figure out whether this face will be stored as f 
            // or ~f, based on whether its nodes produce a normal vector that 
            // points outward (f) or inward (~f). If the nodes n0, n1, n2 are 
            // equivalent to any even permutation of the nodes of the faces 
            // of our tet, the face is stored as f. Otherwise, we use ~f.
            if (nodes_are_even_permutation_of_tet_face(n0, n1, n2, tet))
              mesh->cell_faces[j] = f;
            else
              mesh->cell_faces[j] = ~f;
printf("Cell %d face %d = %d\n", cell, j-mesh->cell_face_offsets[cell], mesh->cell_faces[j]);
            break;
          }
          ++j;
        }
      }
    }
    int_unordered_set_free(intersect01);
    int_unordered_set_free(intersect012);
    int_ptr_unordered_map_free(cells_for_node);
#endif
  }
  mesh->storage->cell_face_capacity = 4*mesh->num_cells;

  // Compute the mesh's geometry.
  mesh_compute_geometry(mesh);

  // Clean up.
  free(nodes);
  free(faces);
  free(tets);
  int_table_free(edge_for_nodes);

  return mesh;
}

