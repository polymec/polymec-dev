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
#include "core/slist.h"
#include "geometry/create_tetgen_mesh.h"

typedef struct
{
  int num_nodes; // 4 for order 1, 10 for order 2.
  int nodes[10];  
  int attribute; // 0 for none, positive for actual attribute.
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

static point_t* read_nodes(FILE* node_file, int* num_nodes)
{
  *num_nodes = -1;
  point_t* nodes = NULL;
  char line[4096];
  int nodes_read = 0;
  while (fscanf(node_file, "%s\n", line))
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
      {
        fclose(node_file);
        polymec_error("Node file has bad header.");
      }
      if (dim != 3)
      {
        fclose(node_file);
        polymec_error("Node file is not 3-dimensional.");
      }
      nodes = malloc(sizeof(point_t) * (*num_nodes));
      continue;
    }

    // Read the coordinates of the next node
    char garbage[4096];
    int node_id;
    int num_items = sscanf(line, "%d %lg %lg %lg%s\n", &node_id, 
                           &nodes[nodes_read].x, &nodes[nodes_read].y, 
                           &nodes[nodes_read].z, garbage);
    if (num_items < 4)
      polymec_error("Bad line in nodes file after %d nodes read.\n", nodes_read);

    ++nodes_read;
    if (nodes_read == *num_nodes) break;
  }
  if (nodes_read != *num_nodes)
  {
    fclose(node_file);
    polymec_error("Node file claims to contains %d nodes, but %d were read.", *num_nodes, nodes_read);
  }
  return nodes;
}

static tet_t* read_tets(FILE* tet_file, int* num_tets)
{
  tet_t* tets;
  return tets;
}

static tet_face_t* read_faces(FILE* face_file, int* num_faces)
{
  tet_face_t* faces;
  return faces;
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

mesh_t* create_tetgen_mesh(MPI_Comm comm, 
                           const char* node_file,
                           const char* ele_file,
                           const char* face_file)
{
  FILE* node_fp = fopen(node_file, "r");
  if (node_fp == NULL)
    polymec_error("TetGen node file '%s' not found.", node_file);
  int num_nodes;
  point_t* nodes = read_nodes(node_fp, &num_nodes);
  fclose(node_fp);

  FILE* ele_fp = fopen(ele_file, "r");
  if (ele_fp == NULL)
    polymec_error("TetGen element file '%s' not found.", ele_file);
  int num_tets;
  tet_t* tets = read_tets(ele_fp, &num_tets);
  fclose(ele_fp);

  FILE* face_fp = fopen(face_file, "r");
  if (face_fp == NULL)
    polymec_error("TetGen face file '%s' not found.", face_file);
  int num_faces;
  tet_face_t* faces = read_faces(face_fp, &num_faces);
  fclose(face_fp);

  // Compute the number of edges, which we aren't given by polytope.
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

#if 0
  // Face <-> node connectivity.
  memcpy(mesh->face_node_offsets, tess->face_offsets, sizeof(int) * (tess->num_faces+1));
  mesh->face_nodes = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int) * tess->face_offsets[tess->num_faces], 0);
  memcpy(mesh->face_nodes, tess->face_nodes, sizeof(int) * tess->face_offsets[tess->num_faces]);
  mesh->storage->face_node_capacity = tess->cell_offsets[tess->num_cells];

  // Face <-> edge connectivity.
  memset(mesh->face_edge_offsets, 0, sizeof(int) * (mesh->num_faces + 1));
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int Ne = tess->face_offsets[f+1] - tess->face_offsets[f];
    mesh->face_edge_offsets[f+1] = Ne;
    for (int e = 0; e < Ne; ++e)
    {
      int offset = tess->face_offsets[f];
      int n1 = (int)tess->face_nodes[offset+e];
      int n2 = (int)tess->face_nodes[offset+(e+1)%Ne];
      int edge_id = *int_table_get(edge_for_nodes, n1, n2);
      mesh->storage->face_edge_capacity = round_to_pow2(edge_id+1);
      mesh->face_edges = ARENA_REALLOC(mesh->arena, mesh->face_edges, sizeof(int) * mesh->storage->face_edge_capacity, 0);
      mesh->face_edges[tess->face_offsets[f+1]+e] = edge_id;
    }
  }

  // Cell <-> face connectivity.
  memcpy(mesh->cell_face_offsets, tess->cell_offsets, sizeof(int) * (tess->num_cells+1));
  mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int) * tess->cell_offsets[tess->num_cells], 0);
  memcpy(mesh->cell_faces, tess->cell_faces, sizeof(int) * tess->cell_offsets[tess->num_cells]);
  memcpy(mesh->face_cells, tess->face_cells, sizeof(int) * 2 * tess->num_faces);
  mesh->storage->cell_face_capacity = tess->cell_offsets[tess->num_cells];
#endif

  // Compute the mesh's geometry.
  mesh_compute_geometry(mesh);

  // Clean up.
  int_table_free(edge_for_nodes);

  return mesh;
}

