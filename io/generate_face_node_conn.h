#ifndef POLYMEC_GENERATE_FACE_NODE_CONN_H
#define POLYMEC_GENERATE_FACE_NODE_CONN_H

#include "core/io.h"

#ifdef __cplusplus
extern "C" {
#endif

// Given a mesh, generates face-node connectivity information.
// face_nodes - A pointer that will be set to a newly allocated array 
//              storing all nodes of faces in traversal order, starting 
//              with the first node of the 0th face and ending with the 
//              last node of the last face. This array must be freed with 
//              free.
// face_node_offsets - An array of length (mesh->num_faces+1) whose ith 
//                     element is the index in face_nodes containing the 
//                     first node of face i.
void generate_face_node_conn(mesh_t* mesh,
                             int** face_nodes,
                             int* face_node_offsets);

#ifdef __cplusplus
}
#endif

#endif
