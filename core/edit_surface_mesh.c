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

int surface_mesh_add_node(surface_mesh_t* surface)
{
  if (surface->num_nodes+1 > surface->storage->node_capacity)
  {
    while (surface->num_nodes+1 > surface->storage->node_capacity)
      surface->storage->node_capacity *= 2;
    surface->nodes = arena_realloc(surface->arena, surface->nodes, sizeof(node_t)*surface->storage->node_capacity, 0);
  }
  surface->num_nodes++;
  return surface->num_nodes-1;
}

void surface_mesh_delete_node(surface_mesh_t* surface, int i)
{
  // Swap the ith node with the end.
  if (i < surface->num_nodes)
  {
    int last = surface->num_nodes - 1;
    node_t tmp = surface->nodes[last];
    surface->nodes[last] = surface->nodes[i];
    surface->nodes[i] = tmp;
  }
}

int surface_mesh_add_edge(surface_mesh_t* surface)
{
  if (surface->num_edges+1 > surface->storage->edge_capacity)
  {
    while (surface->num_edges+1 > surface->storage->edge_capacity)
      surface->storage->edge_capacity *= 2;
    surface->edges = arena_realloc(surface->arena, surface->edges, sizeof(edge_t)*surface->storage->edge_capacity, 0);
  }
  surface->num_edges++;
  return surface->num_edges-1;
}

void surface_mesh_delete_edge(surface_mesh_t* surface, int i)
{
  // Swap the ith edge with the end.
  if (i < surface->num_edges)
  {
    int last = surface->num_edges - 1;
    edge_t tmp = surface->edges[last];
    surface->edges[last] = surface->edges[i];
    surface->edges[i] = tmp;
  }
}

int surface_mesh_add_face(surface_mesh_t* surface)
{
  if (surface->num_faces+1 > surface->storage->face_capacity)
  {
    while (surface->num_faces+1 > surface->storage->face_capacity)
      surface->storage->face_capacity *= 2;
    surface->faces = arena_realloc(surface->arena, surface->faces, sizeof(face_t)*surface->storage->face_capacity, 0);
  }
  surface->num_faces++;
  return surface->num_faces-1;
}

void surface_mesh_delete_face(surface_mesh_t* surface, int i)
{
  // Swap the ith face with the end.
  if (i < surface->num_faces)
  {
    int last = surface->num_faces - 1;
    face_t tmp = surface->faces[last];
    surface->faces[last] = surface->faces[i];
    surface->faces[i] = tmp;
  }
}

void surface_mesh_add_edge_to_face(surface_mesh_t* surface, edge_t* edge, face_t* face)
{
  if (face->edges == NULL)
  {
    face->num_edges = 1;
    face->edges = ARENA_MALLOC(surface->arena, sizeof(edge_t*)*4, 0);
    face->edges[0] = edge;
  }
  else
  {
#ifndef NDEBUG
    // Make sure this edge isn't already attached.
    for (int e = 0; e < face->num_edges; ++e)
    {
      ASSERT(edge != face->edges[e]);
    }
#endif
    face->num_edges++;
    int ne = MAX(round_to_pow2(face->num_edges), 4);
    face->edges = ARENA_REALLOC(surface->arena, face->edges, sizeof(edge_t*)*ne, 0);
    face->edges[face->num_edges-1] = edge;
  }
}

void surface_mesh_remove_edge_from_face(surface_mesh_t* surface, edge_t* edge, face_t* face)
{
  for (int e = 0; e < face->num_edges; ++e)
  {
    if (face->edges[e] == edge)
    {
      face->num_edges--;
      face->edges[e] = face->edges[face->num_edges];
      break;
    }
  }
}

#ifdef __cplusplus
}
#endif

