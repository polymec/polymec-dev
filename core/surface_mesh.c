#include "core/surface_mesh.h"
#include "core/mesh_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

surface_mesh_t* surface_mesh_new(int num_faces, int num_edges, int num_nodes)
{
  ARENA* a = arena_open(&arena_defaults, 0);
  surface_mesh_t* s = surface_mesh_new_with_arena(a, num_faces, num_edges, num_nodes);
  s->close_arena = true;
  return s;
}

surface_mesh_t* surface_mesh_new_with_arena(ARENA* arena, int num_faces, int num_edges, int num_nodes)
{
  ASSERT(num_faces > 0);
  ASSERT(num_edges > 0);
  ASSERT(num_nodes > 0);

  surface_mesh_t* s = ARENA_MALLOC(arena, sizeof(surface_mesh_t), 0);
  s->arena = arena;
  s->close_arena = false;

  // NOTE: We round stored elements up to the nearest power of 2.
  int face_cap = round_to_pow2(num_faces);
  s->faces = ARENA_MALLOC(s->arena, sizeof(face_t)*face_cap, 0);
  s->num_faces = num_faces;

  int edge_cap = round_to_pow2(num_edges);
  s->edges = ARENA_MALLOC(s->arena, sizeof(edge_t)*edge_cap, 0);
  s->num_edges = num_edges;

  int node_cap = round_to_pow2(num_nodes);
  s->nodes = ARENA_MALLOC(s->arena, sizeof(node_t)*node_cap, 0);
  s->num_nodes = num_nodes;

  // Storage information.
  s->storage = mesh_storage_new_with_arena(arena);
  s->storage->node_capacity = node_cap;
  s->storage->edge_capacity = edge_cap;
  s->storage->face_capacity = face_cap;
  s->storage->cell_capacity = 0;

  return s;
}

void surface_mesh_free(surface_mesh_t* surface)
{
  ASSERT(surface != NULL);

  for (int i = 0; i < surface->num_faces; ++i)
  {
    if (surface->faces[i].edges != NULL)
      ARENA_FREE(surface->arena, surface->faces[i].edges);
  }
  ARENA_FREE(surface->arena, surface->faces);
  ARENA_FREE(surface->arena, surface->edges);
  ARENA_FREE(surface->arena, surface->nodes);

  mesh_storage_free(surface->storage);

  ARENA* arena = surface->arena;
  bool close_arena = surface->close_arena;
  ARENA_FREE(arena, surface);
  if (close_arena)
    arena_close(arena);
}

#ifdef __cplusplus
}
#endif

