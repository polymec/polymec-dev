#ifndef POLYMEC_FACETED_SURFACE_H
#define POLYMEC_FACETED_SURFACE_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// This data type represents a faceted surface, consisting of faces 
// connected with edges and nodes in three-dimensional space.
typedef struct 
{
  // Faces, indexed from 0 to F-1.
  face_t* faces;
  // Total number of faces in the surface.
  int num_faces;
  // Edges, indexed from 0 to E-1.
  edge_t* edges;
  // Total number of edges in the surface.
  int num_edges;
  // Nodes, indexed from 0 to N-1.
  node_t* nodes;
  // Total number of nodes in the surface.
  int num_nodes;

  // Storage information -- used internally.
  ARENA* arena;
  bool close_arena;
  mesh_storage_t* storage;
} faceted_surface_t;

// Construct a new surface with the given number of faces, edges, and nodes. 
// This function does not provide any description of the mesh's topology and 
// is only useful in the construction of mesh generation algorithms.
faceted_surface_t* faceted_surface_new(int num_faces, int num_edges, int num_nodes);

// Construct a new surface, using the given arena for memory allocations.
faceted_surface_t* faceted_surface_new_with_arena(ARENA* arena, int num_faces, int num_edges, int num_nodes);

// Destroys the given surface.
void faceted_surface_free(faceted_surface_t* surface);

#ifdef __cplusplus
}
#endif

#endif

