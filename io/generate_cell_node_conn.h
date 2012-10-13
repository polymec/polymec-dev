#ifndef ARBI_GENERATE_CELL_NODE_CONN_H
#define ARBI_GENERATE_CELL_NODE_CONN_H

#include "core/io.h"

#ifdef __cplusplus
extern "C" {
#endif

// Given a mesh and existing face-node connectivity data, generates cell-node 
// connectivity information.
// face_nodes - An array storing all nodes of faces in traversal 
//              order, starting with the first node of the 0th face and ending 
//              with the last node of the last face.
// face_node_offsets - An array of length (mesh->num_faces+1) whose ith 
//                     element is the index in face_nodes containing the 
//                     first node of face i.
// cell_nodes - A pointer that will be set to a newly allocated array 
//              storing all nodes of cells, starting with the first node of 
//              the 0th cell and ending with the last node of the last cell. 
//              This array must be freed with free.
// cell_node_offsets - An array of length (mesh->num_cells+1) whose ith 
//                     element is the index in cell_nodes containing the 
//                     first node of cell i.
void generate_cell_node_conn(mesh_t* mesh,
                             int* face_nodes,
                             int* face_node_offsets,
                             int** cell_nodes,
                             int* cell_node_offsets);

#ifdef __cplusplus
}
#endif

#endif
