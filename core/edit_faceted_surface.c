#include "core/edit_faceted_surface.h"
#include "core/mesh_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

int add_faceted_surface_node(faceted_surface_t* surface)
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

void delete_faceted_surface_node(faceted_surface_t* surface, int i)
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

int add_faceted_surface_edge(faceted_surface_t* surface)
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

void delete_faceted_surface_edge(faceted_surface_t* surface, int i)
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

int add_faceted_surface_face(faceted_surface_t* surface)
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

void delete_faceted_surface_face(faceted_surface_t* surface, int i)
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

#ifdef __cplusplus
}
#endif

