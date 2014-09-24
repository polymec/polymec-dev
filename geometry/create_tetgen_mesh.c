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

#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/text_buffer.h"
#include "core/array_utils.h"
#include "core/partition_mesh.h"
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

static point_t* read_nodes(const char* node_file, int* num_nodes)
{
  *num_nodes = -1;
  point_t* nodes = NULL;

  text_buffer_t* buffer = text_buffer_from_file(node_file);
  if (buffer == NULL)
    polymec_error("TetGen node file '%s' not found.", node_file);
  int nodes_read = 0, pos = 0, line_length;
  char* line;
  while (text_buffer_next(buffer, &pos, &line, &line_length))
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
      if (*num_nodes <= 0)
        polymec_error("Node file has bad number of nodes: %d.", *num_nodes);
      if (dim != 3)
        polymec_error("Node file is not 3-dimensional.");
      nodes = polymec_malloc(sizeof(point_t) * (*num_nodes));
      continue;
    }

    // Read the coordinates of the next node
    char garbage[4096];
    int node_id;
    int num_items = sscanf(line, "%d %lg %lg %lg%s\n", &node_id, 
                           &nodes[nodes_read].x, &nodes[nodes_read].y, 
                           &nodes[nodes_read].z, garbage);
    if (node_id != (nodes_read+1))
      polymec_error("Bad node ID after %d nodes read: %d.", nodes_read, node_id);
    if (num_items < 4)
      polymec_error("Bad line in nodes file after %d nodes read.", nodes_read);
    if (*num_nodes <= 0)
      polymec_error("Bad number of nodes in node file: %d.", *num_nodes);

    ++nodes_read;
    if (nodes_read == *num_nodes) break;
  }
  text_buffer_free(buffer);

  if (nodes_read != *num_nodes)
    polymec_error("Node file claims to contain %d nodes, but %d were read.", *num_nodes, nodes_read);
  return nodes;
}

static tet_t* read_tets(const char* tet_file, int* num_tets)
{
  *num_tets = -1;
  tet_t* tets = NULL;
  text_buffer_t* buffer = text_buffer_from_file(tet_file);
  if (buffer == NULL)
    polymec_error("TetGen element file '%s' not found.", tet_file);
  int tets_read = 0, nodes_per_tet = 4, region_attribute = 0, pos = 0, line_length;
  char* line;
  while (text_buffer_next(buffer, &pos, &line, &line_length))
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
      tets = polymec_malloc(sizeof(tet_t) * (*num_tets));
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
        polymec_error("Bad line in element file after %d tets read.", tets_read);
      if (num_items == 6)
      {
        if (string_is_number(attr_str))
          tet->attribute = atoi(attr_str);
      }
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
        polymec_error("Bad line in element file after %d tets read.", tets_read);
      if (num_items == 12)
      {
        if (string_is_number(attr_str))
          tet->attribute = atoi(attr_str);
      }
    }
    if (tet_id != (tets_read+1))
      polymec_error("Bad tet ID after %d tets read: %d.", tets_read, tet_id);
    ++tets_read;
    if (tets_read == *num_tets) break;
  }
  text_buffer_free(buffer);
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
  text_buffer_t* buffer = text_buffer_from_file(face_file);
  if (buffer == NULL)
    polymec_error("TetGen face file '%s' not found.", face_file);
  int faces_read = 0, boundary_marker = 0, pos = 0, line_length;
  char* line;
  while (text_buffer_next(buffer, &pos, &line, &line_length))
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
      faces = polymec_malloc(sizeof(tet_face_t) * (*num_faces));
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
        polymec_error("Bad line in face file after %d faces read.", faces_read);
      if (num_items == 5)
      {
        if (string_is_number(marker_str))
          face->boundary_marker = atoi(marker_str);
      }
    }
    else
    {
      face->num_nodes = 6;
      int num_items = sscanf(line, "%d %d %d %d %d %d %d%s\n", 
                             &face_id, &face->nodes[0], &face->nodes[1], 
                             &face->nodes[2], &face->nodes[3], &face->nodes[4],
                             &face->nodes[5], marker_str);
      if (num_items < 7)
        polymec_error("Bad line in face file after %d faces read.", faces_read);
      if (num_items == 8)
      {
        if (string_is_number(marker_str))
          face->boundary_marker = atoi(marker_str);
      }
    }
    if (face_id != (faces_read+1))
      polymec_error("Bad face ID after %d faces read: %d.", faces_read, face_id);
    ++faces_read;
    if (faces_read == *num_faces) break;
  }
  text_buffer_free(buffer);
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
  text_buffer_t* buffer = text_buffer_from_file(neigh_file);
  if (buffer == NULL)
    polymec_error("TetGen neighbor file '%s' not found.", neigh_file);
  int tets_read = 0, pos = 0, line_length, num_entries = -1;
  char* line;
  while (text_buffer_next(buffer, &pos, &line, &line_length))
  {
    // Skip lines starting with #.
    if (line[0] == '#') continue;

    // Look for the header if we haven't already read it.
    if (num_entries == -1)
    {
      int four;
      int num_items = sscanf(line, "%d %d\n", &num_entries, &four);
      if (num_items != 2)
        polymec_error("Neighbor file has bad header.");
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
        polymec_error("Bad line in neighbors file after %d tets read.", tets_read);
    if (tet_id != (tets_read+1))
      polymec_error("Bad tet ID after %d tet read: %d.", tets_read, tet_id);
    ++tets_read;
    if (tets_read == num_tets) break;
  }
  text_buffer_free(buffer);
  if (tets_read != num_tets)
    polymec_error("Neighbor file has %d tets, but needs %d.", tets_read, num_tets);

  // TetGen's indices are 1-based, so correct them.
  for (int t = 0; t < num_tets; ++t)
  {
    tet_t* tet = &tets[t];
    for (int n = 0; n < 4; ++n)
    {
      if (tet->neighbors[n] > 0)
        tet->neighbors[n] -= 1;
    }
  }
}

static bool face_points_outward(tet_face_t* face,
                                tet_t* tet,
                                point_t* nodes)
{
  // Compute the face's center.
  point_t* n1 = &nodes[face->nodes[0]];
  point_t* n2 = &nodes[face->nodes[1]];
  point_t* n3 = &nodes[face->nodes[2]];
  point_t xf = {.x = (n1->x + n2->x + n3->x)/3.0,
                .y = (n1->y + n2->y + n3->y)/3.0,
                .z = (n1->z + n2->z + n3->z)/3.0};

  // Compute the face's normal vector, assuming that its nodes are ordered 
  // counterclockwise.
  vector_t x12, x13, nf;
  point_displacement(n1, n2, &x12);
  point_displacement(n1, n3, &x13);
  vector_cross(&x12, &x13, &nf);
  ASSERT(vector_mag(&nf) != 0.0);

  // Get the remaining tetrahedron node.
  point_t* n4 = ((tet->nodes[0] != face->nodes[0]) &&
                 (tet->nodes[0] != face->nodes[1]) && 
                 (tet->nodes[0] != face->nodes[2])) ? &nodes[tet->nodes[0]] :
                ((tet->nodes[1] != face->nodes[0]) &&
                 (tet->nodes[1] != face->nodes[1]) && 
                 (tet->nodes[1] != face->nodes[2])) ? &nodes[tet->nodes[1]] :
                ((tet->nodes[2] != face->nodes[0]) &&
                 (tet->nodes[2] != face->nodes[1]) && 
                 (tet->nodes[2] != face->nodes[2])) ? &nodes[tet->nodes[2]] :
                 &nodes[tet->nodes[3]];

  // Determine whether nf points in the same direction as (xf - n4).
  vector_t d;
  point_displacement(n4, &xf, &d);
  return (vector_dot(&d, &nf) > 0.0);
}

mesh_t* create_tetgen_mesh(MPI_Comm comm,
                           const char* node_file,
                           const char* ele_file,
                           const char* face_file,
                           const char* neigh_file)
{
  int rank, nproc;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nproc);

  mesh_t* mesh = NULL;

  if (rank == 0)
  {
    // Read the information in the file.
    int num_nodes;
    point_t* nodes = read_nodes(node_file, &num_nodes);

    int num_tets;
    tet_t* tets = read_tets(ele_file, &num_tets);

    int num_faces;
    int nodes_per_face = (tets[0].num_nodes == 4) ? 3 : 6;
    tet_face_t* faces = read_faces(face_file, nodes_per_face, &num_faces);

    read_neighbors(neigh_file, tets, num_tets);

    // Create a mesh full of tetrahedra (4 faces per cell, 3 nodes per face).
    mesh = mesh_new_with_cell_type(MPI_COMM_SELF, num_tets, 0, num_faces, num_nodes, 4, nodes_per_face);

    // Copy node coordinates.
    memcpy(mesh->nodes, nodes, sizeof(point_t) * num_nodes);

    // Actual connectivity.
    int_tuple_int_unordered_map_t* face_for_nodes = int_tuple_int_unordered_map_new();
    for (int f = 0; f < num_faces; ++f)
    {
      tet_face_t* face = &faces[f];
      for (int n = 0; n < nodes_per_face; ++n)
        mesh->face_nodes[nodes_per_face*f+n] = face->nodes[n];

      // Associate the 3 "primal" nodes with this face.
      int* primal_nodes = int_tuple_new(3);
      for (int i = 0; i < 3; ++i)
        primal_nodes[i] = face->nodes[i];
      int_qsort(primal_nodes, 3);
      int_tuple_int_unordered_map_insert_with_k_dtor(face_for_nodes, primal_nodes, f, int_tuple_free);
    }

    // Cell <-> face connectivity.
    for (int c = 0; c < mesh->num_cells; ++c)
    {
      for (int f = mesh->cell_face_offsets[c]; f < mesh->cell_face_offsets[c+1]; ++f)
        mesh->cell_faces[f] = -1;
    }
    for (int f = 0; f < num_faces; ++f)
    {
      mesh->face_cells[2*f]   = -1;
      mesh->face_cells[2*f+1] = -1;
    }
    {
      // Use a triple for querying faces.
      int* nodes = int_tuple_new(3);

      // Loop over cells and find the faces connecting them to their neighbors.
      for (int c = 0; c < mesh->num_cells; ++c)
      {
        tet_t* t = &tets[c];

        // Tet face-node table.
        static int tet_face_nodes[4][3] = {{1, 2, 3},  // face 1 has nodes 2, 3, 4
          {2, 0, 3},  // face 2 has nodes 3, 1, 4
          {0, 1, 3},  // face 3 has nodes 1, 2, 4
          {0, 1, 2}}; // face 4 has nodes 1, 2, 3

        // Figure out each of the connections by examining their common nodes.
        // We use TetGen's indexing scheme (see TetGen documentation), which 
        // states that neighbor n of a tet shares the face of that tet that 
        // is opposite of node n in the tet.
        for (int n = 0; n < 4; ++n)
        {
          // Nodes of cell c on this face.
          nodes[0] = t->nodes[tet_face_nodes[n][0]];
          nodes[1] = t->nodes[tet_face_nodes[n][1]];
          nodes[2] = t->nodes[tet_face_nodes[n][2]];
          int_qsort(nodes, 3);

          // Find the face with these nodes.
          int* face_p = int_tuple_int_unordered_map_get(face_for_nodes, nodes);
          if (face_p == NULL)
            polymec_error("TetGen files are inconsistent (cell %d does not have a face with nodes %d, %d, %d)", c+1, nodes[0]+1, nodes[1]+1, nodes[2]+1);
          int face = *face_p;

          // Determine whether the face has an outward or inward normal w.r.t. 
          // the cell.
          tet_face_t* tf = &faces[face];
          bool outward_normal = face_points_outward(tf, t, mesh->nodes);

          // Get the neighbor tet.
          int cn = t->neighbors[n];
          if (cn == -1)
          {
            // Set up the face.
            mesh->cell_faces[mesh->cell_face_offsets[c]+n] = outward_normal ? face : ~face;
            mesh->face_cells[2*face] = c;
          }
          else if (cn > c)
          {
            tet_t* tn = &tets[cn];

            // Find the neighbor index of c within cn.
            int n1 = (tn->neighbors[0] == c) ? 0 :
              (tn->neighbors[1] == c) ? 1 : 
              (tn->neighbors[2] == c) ? 2 : 3;

            // Associate the face with both of these cells.
            mesh->cell_faces[mesh->cell_face_offsets[c]+n]  = outward_normal ? face : ~face;
            mesh->cell_faces[mesh->cell_face_offsets[cn]+n1] = outward_normal ? ~face : face;

            // Associate the cells with the face.
            mesh->face_cells[2*face]   = c;
            mesh->face_cells[2*face+1] = cn;
          }
        }
      }

      // Clean up.
      int_tuple_free(nodes);
    }

    // Build edges.
    mesh_construct_edges(mesh);

    // Compute the mesh's geometry.
    mesh_compute_geometry(mesh);

    // Set up tags for faces and cells.
    static const int max_num_attr = 1024;
    int boundary_markers[max_num_attr], attributes[max_num_attr];
    for (int i = 0; i < max_num_attr; ++i)
      boundary_markers[i] = attributes[i] = 0;
    for (int f = 0; f < num_faces; ++f)
    {
      ASSERT(faces[f].boundary_marker < max_num_attr);
      if (faces[f].boundary_marker != -1)
        boundary_markers[faces[f].boundary_marker]++;
    }
    int* face_tags[max_num_attr];
    for (int i = 0; i < max_num_attr; ++i)
    {
      if (boundary_markers[i] > 0)
      {
        char tag_name[16];
        snprintf(tag_name, 16, "%d", i);
        face_tags[i] = mesh_create_tag(mesh->face_tags, tag_name, boundary_markers[i]);
      }
    }
    memset(boundary_markers, 0, sizeof(int) * max_num_attr);
    for (int f = 0; f < num_faces; ++f)
    {
      int m = faces[f].boundary_marker;
      if (m != -1)
      {
        face_tags[m][boundary_markers[m]] = f;
        boundary_markers[m]++;
      }
    }

    for (int t = 0; t < num_tets; ++t)
    {
      // If this is a "normal" attribute, we assign it to a tag.
      if ((tets[t].attribute < max_num_attr) && (tets[t].attribute != -1))
        attributes[tets[t].attribute]++;
      // Otherwise it's probably something to do with adaptive resolution.
    }
    int* cell_tags[max_num_attr];
    for (int i = 0; i < max_num_attr; ++i)
    {
      if (attributes[i] > 0)
      {
        char tag_name[16];
        snprintf(tag_name, 16, "%d", i);
        cell_tags[i] = mesh_create_tag(mesh->cell_tags, tag_name, attributes[i]);
      }
    }
    memset(attributes, 0, sizeof(int) * max_num_attr);
    for (int t = 0; t < num_tets; ++t)
    {
      int a = tets[t].attribute;
      if ((a < max_num_attr) && (a != -1))
      {
        cell_tags[a][attributes[a]] = t;
        attributes[a]++;
      }
    }

    // Clean up.
    polymec_free(nodes);
    polymec_free(faces);
    polymec_free(tets);
    int_tuple_int_unordered_map_free(face_for_nodes);
  }

  // Partition the mesh (without weights).
  partition_mesh(&mesh, comm, NULL, 0.05);
  
  mesh_add_feature(mesh, TETRAHEDRAL);
  return mesh;
}

