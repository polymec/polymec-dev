#ifndef POLYMEC_MESH_DELTA_H
#define POLYMEC_MESH_DELTA_H

#include "core/mesh.h"

// A mesh_delta is a change to the topology of the mesh. Objects of this 
// type are garbage-collected.
typedef struct mesh_delta_t mesh_delta_t;

// Applies the given delta to the given mesh.
void mesh_delta_apply(mesh_delta_t* delta, mesh_t* mesh);

// Constructs a mesh_delta representing the inverse of this one.
mesh_delta_t* mesh_delta_inverse(mesh_delta_t* delta);

// Constructs a mesh_delta that swaps the positions of two elements of the 
// given type within the corresponding array within the mesh.
mesh_delta_t* swap_mesh_delta_new(mesh_centering_t type, int index1, int index2);

// Constructs a mesh_delta that appends a new element of the given type to 
// the mesh's corresponding array. This works for all types except MESH_NODE,
// which requires a set of coordinates for the added node.
mesh_delta_t* append_mesh_delta_new(mesh_centering_t type);

// Constructs a mesh delta that appends a new node to the mesh with the given 
// coordinates.
mesh_delta_t* append_node_mesh_delta_new(point_t* x);

// Constructs a mesh_delta that pops an element of the given type off the 
// back of the mesh's corresponding array.
mesh_delta_t* pop_mesh_delta_new(mesh_centering_t type);

// Constructs a mesh_delta that attachs an element of the given type to 
// a given parent in the mesh.
mesh_delta_t* attach_mesh_delta_new(mesh_centering_t type, int index, int parent_index);

// Constructs a mesh_delta that detachs an element of the given type from 
// a given parent in the mesh.
mesh_delta_t* detach_mesh_delta_new(mesh_centering_t type, int index, int parent_index);

// Writes a text representation of the mesh_delta to the given file.
void mesh_delta_fprintf(mesh_delta_t* delta, FILE* file);

#endif

