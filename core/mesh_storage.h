#ifndef POLYMEC_MESH_STORAGE_H
#define POLYMEC_MESH_STORAGE_H

#include "core/mesh.h"
#include "arena/proto.h"

struct mesh_storage_t
{
  ARENA* arena;
  bool close_arena;

  int node_capacity;
  int edge_capacity;
  int face_capacity;
  int cell_capacity;
};

// Initializes a new storage mechanism for a mesh.
mesh_storage_t* mesh_storage_new();

// Initializes a new storage mechanism for a mesh with the given ARENA.
mesh_storage_t* mesh_storage_new_with_arena(ARENA* arena);

// Frees the given storage mechanism.
void mesh_storage_free(mesh_storage_t* storage);

#endif

