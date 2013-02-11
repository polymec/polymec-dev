#include "core/surface_mesh.h"
#include "core/edit_surface_mesh.h"
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

surface_mesh_t* surface_mesh_from_bbox(bbox_t* bbox)
{
  ARENA* a = arena_open(&arena_defaults, 0);
  surface_mesh_t* s = surface_mesh_from_bbox_with_arena(a, bbox);
  s->close_arena = true;
  return s;
}

surface_mesh_t* surface_mesh_from_bbox_with_arena(ARENA* arena, bbox_t* bbox)
{
  // The bounding box has six faces, each of which is rectangular. We split 
  // each rectangle into 2 triangles.
  int num_faces = 2*6, num_edges = 12+4, num_nodes = 8;
  surface_mesh_t* s = surface_mesh_new_with_arena(arena, num_faces, num_edges, num_nodes);

  // For this surface mesh, we use the usual finite-element-inspired node 
  // ordering for a hexahedron.

  // Edges.
  s->edges[0].node1 = &s->nodes[0];
  s->edges[0].node2 = &s->nodes[4];

  s->edges[1].node1 = &s->nodes[4];
  s->edges[1].node2 = &s->nodes[7];

  s->edges[2].node1 = &s->nodes[3];
  s->edges[2].node2 = &s->nodes[7];

  s->edges[3].node1 = &s->nodes[0];
  s->edges[3].node2 = &s->nodes[3];

  s->edges[4].node1 = &s->nodes[1];
  s->edges[4].node2 = &s->nodes[5];

  s->edges[5].node1 = &s->nodes[5];
  s->edges[5].node2 = &s->nodes[6];

  s->edges[6].node1 = &s->nodes[2];
  s->edges[6].node2 = &s->nodes[6];

  s->edges[7].node1 = &s->nodes[1];
  s->edges[7].node2 = &s->nodes[2];

  s->edges[8].node1 = &s->nodes[0];
  s->edges[8].node2 = &s->nodes[1];

  s->edges[9].node1 = &s->nodes[4];
  s->edges[9].node2 = &s->nodes[5];

  s->edges[10].node1 = &s->nodes[2];
  s->edges[10].node2 = &s->nodes[3];

  s->edges[11].node1 = &s->nodes[6];
  s->edges[11].node2 = &s->nodes[7];

  s->edges[12].node1 = &s->nodes[0];
  s->edges[12].node2 = &s->nodes[7];

  s->edges[13].node1 = &s->nodes[1];
  s->edges[13].node2 = &s->nodes[4];

  s->edges[14].node1 = &s->nodes[2];
  s->edges[14].node2 = &s->nodes[5];

  s->edges[15].node1 = &s->nodes[3];
  s->edges[15].node2 = &s->nodes[6];

  // Triangles 0 and 1 (-x).
  surface_mesh_add_edge_to_face(s, &s->edges[0], &s->faces[0]);
  surface_mesh_add_edge_to_face(s, &s->edges[1], &s->faces[0]);
  surface_mesh_add_edge_to_face(s, &s->edges[12], &s->faces[0]);

  surface_mesh_add_edge_to_face(s, &s->edges[2], &s->faces[1]);
  surface_mesh_add_edge_to_face(s, &s->edges[3], &s->faces[1]);
  surface_mesh_add_edge_to_face(s, &s->edges[12], &s->faces[1]);

  // Triangles 2 and 3 (+x).
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[2]);
  surface_mesh_add_edge_to_face(s, &s->edges[5], &s->faces[2]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[2]);

  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[3]);
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[3]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[3]);

  // Triangles 4 and 5 (-y). // FIXME
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[4]);
  surface_mesh_add_edge_to_face(s, &s->edges[5], &s->faces[4]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[4]);

  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[5]);
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[5]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[5]);

  // Triangles 6 and 7 (+y). // FIXME
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[6]);
  surface_mesh_add_edge_to_face(s, &s->edges[5], &s->faces[6]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[6]);

  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[7]);
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[7]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[7]);

  // Triangles 8 and 9 (-z). // FIXME
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[8]);
  surface_mesh_add_edge_to_face(s, &s->edges[5], &s->faces[8]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[8]);

  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[9]);
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[9]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[9]);

  // Triangles 10 and 11 (+z). // FIXME
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[10]);
  surface_mesh_add_edge_to_face(s, &s->edges[5], &s->faces[10]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[10]);

  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[11]);
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[11]);
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[11]);

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

