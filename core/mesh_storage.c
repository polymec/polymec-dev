// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
  storage->cell_face_capacity = 0;
  storage->face_edge_capacity = 0;
  storage->face_node_capacity = 0;
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

