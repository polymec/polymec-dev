// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_GENERATE_CELL_NODE_CONN_H
#define POLYMEC_GENERATE_CELL_NODE_CONN_H

#include "core/io.h"

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

#endif
