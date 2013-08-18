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

#ifndef POLYMEC_SURFACE_MESH_H
#define POLYMEC_SURFACE_MESH_H

#include "core/mesh.h"

// This data type represents a surface mesh consisting of triangular faces 
// connected with edges and nodes in three-dimensional space.
typedef struct 
{
  // Triangular Faces, indexed from 0 to F-1.
  face_t* faces;
  // Total number of faces in the surface.
  int num_faces;
  // Edges, indexed from 0 to E-1.
  edge_t* edges;
  // Total number of edges in the surface.
  int num_edges;
  // Nodes, indexed from 0 to N-1.
  node_t* nodes;
  // Total number of nodes in the surface.
  int num_nodes;

  // Storage information -- used internally.
  ARENA* arena;
  bool close_arena;
  mesh_storage_t* storage;
} surface_mesh_t;

// Constructs a new surface with the given number of faces, edges, and nodes. 
// This function does not provide any description of the mesh's topology and 
// is only useful in the construction of mesh generation algorithms.
surface_mesh_t* surface_mesh_new(int num_faces, int num_edges, int num_nodes);

// Constructs a new surface, using the given arena for memory allocations.
surface_mesh_t* surface_mesh_new_with_arena(ARENA* arena, int num_faces, int num_edges, int num_nodes);

// Constructs a new surface that represents the given bounding box.
surface_mesh_t* surface_mesh_from_bbox(bbox_t* bbox);

// Constructs a new surface that represents the given bounding box, using the 
// given arena.
surface_mesh_t* surface_mesh_from_bbox_with_arena(ARENA* arena, bbox_t* bbox);

// Destroys the given surface.
void surface_mesh_free(surface_mesh_t* surface);

#endif

