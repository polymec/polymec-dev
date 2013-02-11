#ifndef POLYMEC_EDIT_SURFACE_MESH_H
#define POLYMEC_EDIT_SURFACE_MESH_H

#include "core/edit_mesh.h"
#include "core/surface_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Appends a new node to the surface, returning the new index.
int surface_mesh_add_node(surface_mesh_t* surface);

// Deletes the ith node from the surface. This function simply removes the 
// node from the surface's list of nodes and makes no attempt to provide 
// topological consistency.
void surface_mesh_delete_node(surface_mesh_t* surface, int i);

// Appends a new edge to the surface, returning the new index.
int surface_mesh_add_edge(surface_mesh_t* surface);

// Deletes the ith edge from the surface. This function simply removes the 
// node from the surface's list of nodes and makes no attempt to provide 
// topological consistency.
void surface_mesh_delete_edge(surface_mesh_t* surface, int i);

// Appends a new face to the surface, returning the new index.
int surface_mesh_add_face(surface_mesh_t* surface);

// Deletes the ith face from the surface. This function simply removes the 
// node from the surface's list of nodes and makes no attempt to provide 
// topological consistency.
void surface_mesh_delete_face(surface_mesh_t* surface, int i);

#ifdef __cplusplus
}
#endif

#endif

