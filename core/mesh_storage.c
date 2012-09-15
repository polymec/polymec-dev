#include "core/mesh_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_storage_t* mesh_storage_new()
{
  ARENA* a = arena_open(&arena_defaults, ARENA_STDLIB);
  mesh_storage_t* s = mesh_storage_new_with_arena(a);
  s->close_arena = true;
  return s;
}

mesh_storage_t* mesh_storage_new_with_arena(ARENA* arena)
{
  mesh_storage_t* storage = arena_malloc(arena, sizeof(mesh_storage_t), 0);
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
  arena_free(a, storage);
  if (close_arena)
    arena_close(a);
}

#ifdef __cplusplus
}
#endif

