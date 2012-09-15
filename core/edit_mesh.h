#ifndef ARBI_EDIT_MESH_H
#define ARBI_EDIT_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Appends a new node to the mesh, returning the new index.
int mesh_add_node(mesh_t* mesh);

// Deletes the ith node from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void mesh_delete_node(mesh_t* mesh, int i);

// Appends new edges to the mesh, returning the new index.
int mesh_add_edge(mesh_t* mesh);

// Deletes the ith edge from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void mesh_delete_edge(mesh_t* mesh, int i);

// Appends a new face to the mesh, returning the new index.
int mesh_add_face(mesh_t* mesh);

// Deletes the ith face from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void mesh_delete_face(mesh_t* mesh, int i);

// Appends a new cell to the mesh, returning the new index.
int mesh_add_cell(mesh_t* mesh);

// Deletes the ith cell from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void mesh_delete_cell(mesh_t* mesh, int i);

// Attaches the given edge to the given face.
void mesh_add_edge_to_face(mesh_t* mesh, edge_t* edge, face_t* face);

// Attaches the given face to the given cell.
void mesh_add_face_to_cell(mesh_t* mesh, face_t* face, cell_t* cell);

// Removes the given edge from the given face, or does nothing if the 
// edge is not attached to the face.
void mesh_remove_edge_from_face(mesh_t* mesh, edge_t* edge, face_t* face);

// Removes the given face from the given cell, or does nothing if the 
// face is not attached to the cell.
void mesh_remove_face_from_cell(mesh_t* mesh, face_t* face, cell_t* cell);

#ifdef __cplusplus
}
#endif

#endif

