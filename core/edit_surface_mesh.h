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

#ifndef POLYMEC_EDIT_SURFACE_MESH_H
#define POLYMEC_EDIT_SURFACE_MESH_H

#include "core/edit_mesh.h"
#include "core/surface_mesh.h"

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

// Attaches the given edge to the given face.
void surface_mesh_add_edge_to_face(surface_mesh_t* surface, edge_t* edge, face_t* face);

// Removes the given edge from the given face, or does nothing if the 
// edge is not attached to the face.
void surface_mesh_remove_edge_from_face(surface_mesh_t* surface, edge_t* edge, face_t* face);

#endif

