#ifndef ARBI_EDIT_FACETED_SURFACE_H
#define ARBI_EDIT_FACETED_SURFACE_H

#include "core/edit_mesh.h"
#include "core/faceted_surface.h"

#ifdef __cplusplus
extern "C" {
#endif

// Appends a new node to the surface, returning the new index.
int add_faceted_surface_node(faceted_surface_t* surface);

// Deletes the ith node from the surface. This function simply removes the 
// node from the surface's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_faceted_surface_node(faceted_surface_t* surface, int i);

// Appends a new edge to the surface, returning the new index.
int add_faceted_surface_edge(faceted_surface_t* surface);

// Deletes the ith edge from the surface. This function simply removes the 
// node from the surface's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_faceted_surface_edge(faceted_surface_t* surface, int i);

// Appends a new face to the surface, returning the new index.
int add_faceted_surface_face(faceted_surface_t* surface);

// Deletes the ith face from the surface. This function simply removes the 
// node from the surface's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_faceted_surface_face(faceted_surface_t* surface, int i);

#ifdef __cplusplus
}
#endif

#endif

