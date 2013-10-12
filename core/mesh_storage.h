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

#ifndef POLYMEC_MESH_STORAGE_H
#define POLYMEC_MESH_STORAGE_H

#include "core/mesh.h"
#include "arena/proto.h"

struct mesh_storage_t
{
  ARENA* arena;
  bool close_arena;

  int cell_face_capacity;
  int face_edge_capacity;
  int face_node_capacity;
};

// Initializes a new storage mechanism for a mesh.
mesh_storage_t* mesh_storage_new();

// Initializes a new storage mechanism for a mesh with the given ARENA.
mesh_storage_t* mesh_storage_new_with_arena(ARENA* arena);

// Frees the given storage mechanism.
void mesh_storage_free(mesh_storage_t* storage);

#endif

