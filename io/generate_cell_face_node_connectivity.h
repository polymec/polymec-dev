
#ifndef ARBI_GENERATE_CELL_FACE_NODE_CONNECTIVITY_H
#define ARBI_GENERATE_CELL_FACE_NODE_CONNECTIVITY_H

#include "core/io.h"

#ifdef __cplusplus
extern "C" {
#endif

// Given a mesh, populates the following areas with cell-face-node 
// connectivity information.
// face_node_counts - An array of length mesh->num_faces whose ith 
//                    element is the number of nodes in the ith face.
// all_face_nodes - A pointer that will be set to a newly allocated array 
//                  storing all nodes of faces in traversal order, starting 
//                  with the first node of the 0th face and ending with the 
//                  last node of the last face. This array must be freed with 
//                  free.
// cell_face_counts - An array of length mesh->num_cells whose ith 
//                    element is the number of faces in the ith cell.
// all_face_nodes - A pointer that will be set to a newly allocated array 
//                  storing all faces of cells in traversal order, starting 
//                  with the first face of the 0th cell and ending with the 
//                  last face of the last cell. This array must be freed with 
//                  free.
void generate_cell_face_node_connectivity(mesh_t* mesh,
                                          int* face_node_counts,
                                          int** all_face_nodes,
                                          int* cell_face_counts,
                                          int** all_cell_faces);

#ifdef __cplusplus
}
#endif

#endif
