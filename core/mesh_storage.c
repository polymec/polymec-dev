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

#include "core/mesh_storage.h"

mesh_storage_t* mesh_storage_new()
{
  ARENA* a = arena_open(&arena_defaults, ARENA_STDLIB);
  mesh_storage_t* s = mesh_storage_new_with_arena(a);
  s->close_arena = true;
  return s;
}

mesh_storage_t* mesh_storage_new_with_arena(ARENA* arena)
{
  mesh_storage_t* storage = ARENA_MALLOC(arena, sizeof(mesh_storage_t), 0);
  storage->arena = arena;
  storage->node_capacity = 0;
  storage->edge_capacity = 0;
  storage->face_capacity = 0;
  storage->cell_capacity = 0;
  storage->close_arena = false;
  return storage;
}

void mesh_storage_free(mesh_storage_t* storage)
{
  ARENA* a = storage->arena;
  bool close_arena = storage->close_arena;
  ARENA_FREE(a, storage);
  if (close_arena)
    arena_close(a);
}

