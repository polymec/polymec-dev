#include "core/surface_mesh.h"
#include "core/edit_surface_mesh.h"
#include "core/mesh_storage.h"

// This function rounds the given number up to the nearest power of 2.
static inline int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

// This computes the center of the triangle with given node positions.
static inline void compute_tri_center(node_t* n1, node_t* n2, node_t* n3, point_t* center)
{
  center->x = (n1->x + n2->x + n3->x) / 3.0;
  center->y = (n1->y + n2->y + n3->y) / 3.0;
  center->z = (n1->z + n2->z + n3->z) / 3.0;
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
  ASSERT(bbox->x2 > bbox->x1);
  ASSERT(bbox->y2 > bbox->y1);
  ASSERT(bbox->z2 > bbox->z1);

  // The bounding box has six faces, each of which is rectangular. We split 
  // each rectangle into 2 triangles.
  int num_faces = 2*6, num_edges = 12+4, num_nodes = 8;
  surface_mesh_t* s = surface_mesh_new_with_arena(arena, num_faces, num_edges, num_nodes);

  // For this surface mesh, we use the usual finite-element-inspired node 
  // ordering for a hexahedron.

  // Node positions.
  s->nodes[0].x = bbox->x1;
  s->nodes[0].y = bbox->y1;
  s->nodes[0].z = bbox->z1;

  s->nodes[1].x = bbox->x2;
  s->nodes[1].y = bbox->y1;
  s->nodes[1].z = bbox->z1;

  s->nodes[2].x = bbox->x2;
  s->nodes[2].y = bbox->y2;
  s->nodes[2].z = bbox->z1;

  s->nodes[3].x = bbox->x1;
  s->nodes[3].y = bbox->y2;
  s->nodes[3].z = bbox->z1;

  s->nodes[4].x = bbox->x1;
  s->nodes[4].y = bbox->y1;
  s->nodes[4].z = bbox->z2;

  s->nodes[5].x = bbox->x2;
  s->nodes[5].y = bbox->y1;
  s->nodes[5].z = bbox->z2;

  s->nodes[6].x = bbox->x2;
  s->nodes[6].y = bbox->y2;
  s->nodes[6].z = bbox->z2;

  s->nodes[7].x = bbox->x1;
  s->nodes[7].y = bbox->y2;
  s->nodes[7].z = bbox->z2;

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

  s->edges[16].node1 = &s->nodes[1];
  s->edges[16].node2 = &s->nodes[3];

  s->edges[17].node1 = &s->nodes[5];
  s->edges[17].node2 = &s->nodes[7];

  // Triangles 0 and 1 (-x).
  double x_face_area = (bbox->y2 - bbox->y1) * (bbox->z2 - bbox->z1);
  surface_mesh_add_edge_to_face(s, &s->edges[0], &s->faces[0]);  // 0->4
  surface_mesh_add_edge_to_face(s, &s->edges[1], &s->faces[0]);  // 4->7
  surface_mesh_add_edge_to_face(s, &s->edges[12], &s->faces[0]); // 7->0
  s->faces[0].area = 0.5 * x_face_area;
  compute_tri_center(&s->nodes[0], &s->nodes[4], &s->nodes[7], &s->faces[0].center);

  surface_mesh_add_edge_to_face(s, &s->edges[2], &s->faces[1]);  // 7->3
  surface_mesh_add_edge_to_face(s, &s->edges[3], &s->faces[1]);  // 3->0
  surface_mesh_add_edge_to_face(s, &s->edges[12], &s->faces[1]); // 0->7
  s->faces[1].area = 0.5 * x_face_area;
  compute_tri_center(&s->nodes[0], &s->nodes[3], &s->nodes[7], &s->faces[1].center);

  // Triangles 2 and 3 (+x).
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[2]);  // 2->6
  surface_mesh_add_edge_to_face(s, &s->edges[5], &s->faces[2]);  // 6->5
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[2]); // 5->2
  s->faces[2].area = 0.5 * x_face_area;
  compute_tri_center(&s->nodes[2], &s->nodes[5], &s->nodes[6], &s->faces[2].center);

  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[3]);  // 2->1
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[3]);  // 1->5
  surface_mesh_add_edge_to_face(s, &s->edges[14], &s->faces[3]); // 5->2
  s->faces[3].area = 0.5 * x_face_area;
  compute_tri_center(&s->nodes[1], &s->nodes[2], &s->nodes[5], &s->faces[3].center);

  // Triangles 4 and 5 (-y). 
  double y_face_area = (bbox->x2 - bbox->x1) * (bbox->z2 - bbox->z1);
  surface_mesh_add_edge_to_face(s, &s->edges[4], &s->faces[4]);  // 1->5
  surface_mesh_add_edge_to_face(s, &s->edges[9], &s->faces[4]);  // 5->4
  surface_mesh_add_edge_to_face(s, &s->edges[13], &s->faces[4]); // 4->1
  s->faces[4].area = 0.5 * y_face_area;
  compute_tri_center(&s->nodes[1], &s->nodes[4], &s->nodes[5], &s->faces[4].center);

  surface_mesh_add_edge_to_face(s, &s->edges[8], &s->faces[5]);  // 0->1
  surface_mesh_add_edge_to_face(s, &s->edges[13], &s->faces[5]); // 1->4
  surface_mesh_add_edge_to_face(s, &s->edges[0], &s->faces[5]);  // 4->0
  s->faces[5].area = 0.5 * y_face_area;
  compute_tri_center(&s->nodes[1], &s->nodes[4], &s->nodes[0], &s->faces[5].center);

  // Triangles 6 and 7 (+y). 
  surface_mesh_add_edge_to_face(s, &s->edges[2], &s->faces[6]);  // 3->7
  surface_mesh_add_edge_to_face(s, &s->edges[11], &s->faces[6]); // 7->6
  surface_mesh_add_edge_to_face(s, &s->edges[15], &s->faces[6]); // 6->3
  s->faces[6].area = 0.5 * y_face_area;
  compute_tri_center(&s->nodes[3], &s->nodes[6], &s->nodes[7], &s->faces[6].center);

  surface_mesh_add_edge_to_face(s, &s->edges[10], &s->faces[7]); // 2->3
  surface_mesh_add_edge_to_face(s, &s->edges[15], &s->faces[7]); // 3->6
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[7]);  // 6->2
  s->faces[7].area = 0.5 * y_face_area;
  compute_tri_center(&s->nodes[3], &s->nodes[6], &s->nodes[2], &s->faces[7].center);

  // Triangles 8 and 9 (-z). 
  double z_face_area = (bbox->x2 - bbox->x1) * (bbox->y2 - bbox->y1);
  surface_mesh_add_edge_to_face(s, &s->edges[7], &s->faces[8]);  // 1->2
  surface_mesh_add_edge_to_face(s, &s->edges[10], &s->faces[8]); // 2->3
  surface_mesh_add_edge_to_face(s, &s->edges[16], &s->faces[8]); // 3->1
  s->faces[8].area = 0.5 * z_face_area;
  compute_tri_center(&s->nodes[1], &s->nodes[2], &s->nodes[3], &s->faces[8].center);

  surface_mesh_add_edge_to_face(s, &s->edges[8], &s->faces[9]);  // 0->1
  surface_mesh_add_edge_to_face(s, &s->edges[16], &s->faces[9]); // 1->3
  surface_mesh_add_edge_to_face(s, &s->edges[3], &s->faces[9]);  // 3->0
  s->faces[9].area = 0.5 * z_face_area;
  compute_tri_center(&s->nodes[1], &s->nodes[0], &s->nodes[3], &s->faces[9].center);

  // Triangles 10 and 11 (+z). 
  surface_mesh_add_edge_to_face(s, &s->edges[6], &s->faces[10]);  // 5->6
  surface_mesh_add_edge_to_face(s, &s->edges[11], &s->faces[10]); // 6->7
  surface_mesh_add_edge_to_face(s, &s->edges[17], &s->faces[10]); // 7->5
  s->faces[10].area = 0.5 * z_face_area;
  compute_tri_center(&s->nodes[5], &s->nodes[6], &s->nodes[7], &s->faces[10].center);

  surface_mesh_add_edge_to_face(s, &s->edges[9], &s->faces[11]);  // 4->5
  surface_mesh_add_edge_to_face(s, &s->edges[17], &s->faces[11]); // 5->7
  surface_mesh_add_edge_to_face(s, &s->edges[1], &s->faces[11]);  // 7->4
  s->faces[11].area = 0.5 * z_face_area;
  compute_tri_center(&s->nodes[4], &s->nodes[6], &s->nodes[7], &s->faces[11].center);

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

