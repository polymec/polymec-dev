#ifndef ARBI_EDIT_MESH_H
#define ARBI_EDIT_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Appends a new node to the mesh, returning the new index.
int add_mesh_node(mesh_t* mesh);

// Deletes the ith node from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_mesh_node(mesh_t* mesh, int i);

// Appends a new edge to the mesh, returning the new index.
int add_mesh_edge(mesh_t* mesh);

// Deletes the ith edge from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_mesh_edge(mesh_t* mesh, int i);

// Appends a new face to the mesh, returning the new index.
int add_mesh_face(mesh_t* mesh);

// Deletes the ith face from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_mesh_face(mesh_t* mesh, int i);

// Appends a new cell to the mesh, returning the new index.
int add_mesh_cell(mesh_t* mesh);

// Deletes the ith cell from the mesh. This function simply removes the 
// node from the mesh's list of nodes and makes no attempt to provide 
// topological consistency.
void delete_mesh_cell(mesh_t* mesh, int i);

// Attaches the given edge to the given face.
void add_face_edge(face_t* face, edge_t* edge);

// Attaches the given face to the given cell.
void add_cell_face(cell_t* cell, face_t* face);

// Removes the given edge from the given face, or does nothing if the 
// edge is not attached to the face.
void remove_face_edge(face_t* face, edge_t* edge);

// Removes the given face from the given cell, or does nothing if the 
// face is not attached to the cell.
void remove_cell_face(cell_t* cell, face_t* face);

#ifdef __cplusplus
}
#endif

#endif

