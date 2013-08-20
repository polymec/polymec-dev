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

#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/slist.h"
#include "core/edit_mesh.h"
#include "core/kd_tree.h"
#include "geometry/create_voronoi_mesh.h"
#include "geometry/voronoi_tessellator.h"
#include "geometry/giftwrap_hull.h"

mesh_t* create_voronoi_mesh(point_t* generators, int num_generators, 
                            point_t* ghost_generators, int num_ghost_generators,
                            int_slist_t* deleted_cells)
{
  ASSERT(generators != NULL);
  ASSERT(num_generators >= 2);
  ASSERT(num_ghost_generators >= 0);

  if (deleted_cells != NULL)
    int_slist_clear(deleted_cells);

  // Gather the points to be tessellated.
  int num_points = num_generators + num_ghost_generators;
  point_t* points = malloc(sizeof(point_t) * num_points);
  memcpy(points, generators, sizeof(point_t) * num_generators);
  memcpy(&points[num_generators], ghost_generators, sizeof(point_t) * num_ghost_generators);

  // Perform the tessellation.
  voronoi_tessellator_t* tessellator = voronoi_tessellator_new();
  voronoi_tessellation_t* tessellation = voronoi_tessellator_tessellate(tessellator, points, num_points, deleted_cells);

  // The number of cells generated may be less than the number of generators
  // in cases in which generators lie on planes of the convex hull.
  ASSERT(tessellation->num_cells <= (num_generators + num_ghost_generators));
  free(points);

  // Construct the Voronoi graph.
  mesh_t* mesh = mesh_new(tessellation->num_cells,
                          num_ghost_generators, // ???
                          tessellation->num_faces,
                          tessellation->num_edges,
                          tessellation->num_nodes);
  
  // Node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i].x = tessellation->nodes[3*i];
    mesh->nodes[i].y = tessellation->nodes[3*i+1];
    mesh->nodes[i].z = tessellation->nodes[3*i+2];
  }

  // Edge <-> node connectivity.
  for (int i = 0; i < mesh->num_edges; ++i)
  {
    mesh->edges[i].node1 = &mesh->nodes[tessellation->edges[i].node1];
    int n2 = tessellation->edges[i].node2; 
    ASSERT(n2 >= 0);
    mesh->edges[i].node2 = &mesh->nodes[n2];
  }

  // Face <-> edge connectivity.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    int Ne = tessellation->faces[f].num_edges;
    for (int e = 0; e < Ne; ++e)
    {
      int edge_id = tessellation->faces[f].edges[e];
      mesh_attach_edge_to_face(mesh, &mesh->edges[edge_id], &mesh->faces[f]);
    }
    ASSERT(mesh->faces[f].num_edges == Ne);
  }

  // Cell <-> face connectivity.
  int_ptr_unordered_map_t* oce = int_ptr_unordered_map_new();
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int Nf = tessellation->cells[i].num_faces;
    for (int f = 0; f < Nf; ++f)
    {
      int face_index = tessellation->cells[i].faces[f];
      face_t* face = &mesh->faces[face_index];
      mesh_attach_face_to_cell(mesh, face, &mesh->cells[i]);
      for (int e = 0; e < face->num_edges; ++e)
      {
        int edge_index = tessellation->faces[face_index].edges[e];
        ASSERT((face->edges[e] - &mesh->edges[0]) == edge_index);
      }
    }
  }

  // Compute the mesh's geometry.
  mesh_compute_geometry(mesh);

  // Clean up.
  voronoi_tessellation_free(tessellation);

  // Stick the generators into a point set (kd-tree) that the mesh can 
  // carry with it.
  kd_tree_t* generator_set = kd_tree_new(generators, num_generators);
  mesh_set_property(mesh, "generators", generator_set, DTOR(kd_tree_free));

  return mesh;
}

