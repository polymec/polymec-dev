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

// This implementation of the Voronoi tessellator uses Tetgen. It is not built 
// unless a Tetgen tarball was found in the 3rdparty/ directory.
// Please see the license file in the Tetgen tarball for license information 
// on Tetgen.

#ifndef HAVE_TETGEN
#error "create_voronoi_mesh.c should not be built without Tetgen!"
#endif

#undef HAVE_HDF5 // FIXME (HACK): Avoids warnings with polytope 
#undef HAVE_TETGEN // FIXME (HACK): Avoids warnings with polytope 
#include "polytope_c.h"

#include "core/mesh_storage.h"
#include "core/table.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/slist.h"
#include "core/kd_tree.h"
#include "geometry/create_voronoi_mesh.h"

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

static int_table_t* gather_edges(polytope_tessellation_t* tess, 
                                 int* num_edges)
{
  int_table_t* edge_for_nodes = int_table_new();
  *num_edges = 0;
  for (int f = 0; f < tess->num_faces; ++f)
  {
    int Ne = tess->face_offsets[f+1] - tess->face_offsets[f];
    for (int e = 0; e < Ne; ++e)
    {
      int offset = tess->face_offsets[f];
      int n1 = (int)tess->face_nodes[offset+e];
      int n2 = (int)tess->face_nodes[offset+(e+1)%Ne];
      int_table_insert(edge_for_nodes, n1, n2, *num_edges);
      *num_edges += 1;
    }
  }
  return edge_for_nodes;
}

static mesh_t* mesh_from_unbounded_tessellation(polytope_tessellation_t* tess, 
                                                point_t* generators)
{
  // Compute the number of edges, which we aren't given by polytope.
  int num_edges = 0;
  int_table_t* edge_for_nodes = gather_edges(tess, &num_edges);

  // Create the mesh.
  mesh_t* mesh = mesh_new(tess->num_cells, 0, // ???
                          tess->num_faces,
                          num_edges,
                          tess->num_nodes);
  
  // Copy node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i].x = tess->nodes[3*i];
    mesh->nodes[i].y = tess->nodes[3*i+1];
    mesh->nodes[i].z = tess->nodes[3*i+2];
  }

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

  // Compute the mesh's geometry.
  mesh_compute_geometry(mesh);

  // Clean up.
  int_table_free(edge_for_nodes);

  return mesh;
}

mesh_t* create_voronoi_mesh(point_t* generators, int num_generators, 
                            point_t* ghost_generators, int num_ghost_generators)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 2);
  ASSERT(num_ghost_generators >= 0);

  // Gather the points to be tessellated.
  int num_points = num_generators + num_ghost_generators;
  double* points = malloc(sizeof(point_t) * 3 * num_points);
  for (int i = 0; i < num_generators; ++i)
  {
    points[3*i] = generators[i].x;
    points[3*i+1] = generators[i].y;
    points[3*i+2] = generators[i].z;
  }
  for (int i = 0; i < num_ghost_generators; ++i)
  {
    points[3*(num_generators+i)] = ghost_generators[i].x;
    points[3*(num_generators+i)+1] = ghost_generators[i].y;
    points[3*(num_generators+i)+2] = ghost_generators[i].z;
  }

  // Perform an unbounded tessellation using polytope.
  polytope_tessellator_t* tessellator = tetgen_tessellator_new();
  polytope_tessellation_t* tess = polytope_tessellation_new(3);
  polytope_tessellator_tessellate_unbounded(tessellator, points, num_points, tess);

  // Create a Voronoi mesh from this tessellation, deleting unbounded cells.
  // FIXME: This doesn't accommodate ghost generators.
  mesh_t* mesh = mesh_from_unbounded_tessellation(tess, generators);

  // Clean up.
  polytope_tessellation_free(tess);
  polytope_tessellator_free(tessellator);
  free(points);

  // Stick the generators into a point set (kd-tree) that the mesh can 
  // carry with it.
  kd_tree_t* generator_set = kd_tree_new(generators, num_generators);
  mesh_set_property(mesh, "generators", generator_set, DTOR(kd_tree_free));

  return mesh;
}

static mesh_t* mesh_from_bounded_tessellation(polytope_tessellation_t* tess, 
                                              point_t* generators)
{
  // Compute the number of edges, which we aren't given by polytope.
  int num_edges = 0;
  int_table_t* edge_for_nodes = gather_edges(tess, &num_edges);

  // Create the mesh.
  mesh_t* mesh = mesh_new(tess->num_cells, 0, // ???
                          tess->num_faces,
                          num_edges,
                          tess->num_nodes);
  
  // Copy node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i].x = tess->nodes[3*i];
    mesh->nodes[i].y = tess->nodes[3*i+1];
    mesh->nodes[i].z = tess->nodes[3*i+2];
  }

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
  memcpy(mesh->face_node_offsets, tess->face_offsets, sizeof(int) * (tess->num_faces+1));
  mesh->face_nodes = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int) * tess->face_offsets[tess->num_faces], 0);
  memcpy(mesh->face_nodes, tess->face_nodes, sizeof(int) * tess->face_offsets[tess->num_faces]);

  // Face <-> edge connectivity.
  // FIXME: Need to make sure that mesh->face_edges is big enough to 
  // FIXME: hold all the edges!
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
      mesh->face_edges[tess->face_offsets[f+1]+e] = edge_id;
    }
  }

  // Cell <-> face connectivity.
  memcpy(mesh->cell_face_offsets, tess->cell_offsets, sizeof(int) * (tess->num_cells+1));
  mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int) * tess->cell_offsets[tess->num_cells], 0);
  memcpy(mesh->cell_faces, tess->cell_faces, sizeof(int) * tess->cell_offsets[tess->num_cells]);
  memcpy(mesh->face_cells, tess->face_cells, sizeof(int) * 2 * tess->num_faces);

  // Compute the mesh's geometry.
  mesh_compute_geometry(mesh);

  // Clean up.
  int_table_free(edge_for_nodes);

  return mesh;
}

mesh_t* create_voronoi_mesh_in_box(point_t* generators, int num_generators, 
                                   point_t* ghost_generators, int num_ghost_generators,
                                   bbox_t* bounding_box)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 2);
  ASSERT(num_ghost_generators >= 0);

  // Gather the points to be tessellated.
  int num_points = num_generators + num_ghost_generators;
  double* points = malloc(sizeof(point_t) * 3 * num_points);
  for (int i = 0; i < num_generators; ++i)
  {
    points[3*i] = generators[i].x;
    points[3*i+1] = generators[i].y;
    points[3*i+2] = generators[i].z;
  }
  for (int i = 0; i < num_ghost_generators; ++i)
  {
    points[3*(num_generators+i)] = ghost_generators[i].x;
    points[3*(num_generators+i)+1] = ghost_generators[i].y;
    points[3*(num_generators+i)+2] = ghost_generators[i].z;
  }

  // Perform an unbounded tessellation using polytope.
  polytope_tessellator_t* tessellator = tetgen_tessellator_new();
  polytope_tessellation_t* tess = polytope_tessellation_new(3);
  double low[3] = {bounding_box->x1, bounding_box->y1, bounding_box->z1};
  double high[3] = {bounding_box->x2, bounding_box->y2, bounding_box->z2};
  polytope_tessellator_tessellate_in_box(tessellator, points, num_points, 
                                         low, high, tess);

  // Create a Voronoi mesh from this tessellation.
  // FIXME: This doesn't accommodate ghost generators.
  mesh_t* mesh = mesh_from_bounded_tessellation(tess, generators);

  // Clean up.
  polytope_tessellation_free(tess);
  polytope_tessellator_free(tessellator);
  free(points);

  // Stick the generators into a point set (kd-tree) that the mesh can 
  // carry with it.
  kd_tree_t* generator_set = kd_tree_new(generators, num_generators);
  mesh_set_property(mesh, "generators", generator_set, DTOR(kd_tree_free));

  return mesh;
}

