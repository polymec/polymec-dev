// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_MESH_H
#define POLYMEC_MESH_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/adj_graph.h"
#include "core/serializer.h"
#include "arena/proto.h"

// Mesh centerings.
typedef enum
{
  MESH_NODE,
  MESH_EDGE,
  MESH_FACE,
  MESH_CELL
} mesh_centering_t;

// A tagging mechanism for tagging mesh nodes/edges/faces/cells 
// with attributes. 
#ifndef TAGGER_T
#define TAGGER_T
typedef struct tagger_t tagger_t;
#endif

// This mesh feature indicates that a mesh is tetrahedral.
const char* TETRAHEDRAL;

typedef struct mesh_storage_t mesh_storage_t;

// This data type represents an unstructured mesh, consisting of 
// stationary cells and the faces connecting them.
typedef struct 
{
  // MPI communicator.
  MPI_Comm comm;

  // The total number of (locally-owned) cells in the mesh.
  int num_cells;
  // The number of ghost cells on the local domain in the mesh.
  int num_ghost_cells;
  // The offsets of the sets of faces attached to cells, stored in CRS format.
  int* cell_face_offsets;
  // The indices of faces attached to cells, stored in CRS format.
  int* cell_faces;

  // The total number of local faces in the mesh.
  int num_faces;
  // The offsets of the sets of nodes attached to faces, stored in CRS format.
  int* face_node_offsets;
  // The indices of nodes attached to faces, stored in CRS format.
  int* face_nodes;

  // The offsets of the sets of edges attached to faces, stored in CRS format.
  int* face_edge_offsets;
  // The indices of edges attached to faces, stored in CRS format.
  int* face_edges;

  // The cells attached to the faces in the mesh. Each face has 2 cells, 
  // so the first cell of the ith face is face_cells[2*i] and the second 
  // is face_cells[2*i+1].
  int* face_cells;

  // The total number of local edges in the mesh.
  int num_edges;

  // The nodes of the edges in the mesh. Each edge has 2 nodes, so the 
  // the first node of the ith edge is edge_nodes[2*i] and the second is 
  // edge_nodes[2*i+1].
  int* edge_nodes;

  // The total number of local nodes in the mesh.
  int num_nodes;
  // Coordinates of the mesh nodes, indexed from 0 to N-1.
  point_t* nodes;

  // Geometry information.
  real_t* cell_volumes;
  point_t* cell_centers;
  point_t* face_centers;
  real_t* face_areas;
  vector_t* face_normals;

  // Mesh tagging mechanisms.
  tagger_t* cell_tags;
  tagger_t* face_tags;
  tagger_t* edge_tags;
  tagger_t* node_tags;

  // Mesh storage information -- used internally.
  ARENA* arena;
  bool close_arena;
  mesh_storage_t* storage;
} mesh_t;

// Construct a new mesh with the given number of cells, ghost cells, 
// faces, edges, and nodes. This function does not provide any description
// of the mesh's topology and is only useful in the construction of mesh 
// generation algorithms.
mesh_t* mesh_new(MPI_Comm comm, int num_cells, int num_ghost_cells, 
                 int num_faces, int num_edges, int num_nodes);

// Construct a new mesh, using the given arena for memory allocations.
mesh_t* mesh_new_with_arena(ARENA* arena, MPI_Comm comm, int num_cells, 
                            int num_ghost_cells, int num_faces,
                            int num_edges, int num_nodes);

// Destroys the given mesh.
void mesh_free(mesh_t* mesh);

// Validates the mesh, throwing an error if it is topologically invalid.
void mesh_verify(mesh_t* mesh);

// Associates a named piece of metadata (a "property") with the mesh itself.
// This can be used to store information about (for example) how the mesh 
// was generated, which can sometimes be useful. A destructor function can be 
// passed in to handle freeing of resources. If the given property exists 
// on the mesh, it is overwritten.
void mesh_set_property(mesh_t* mesh, const char* property, void* data, void (*dtor)(void*));

// Retrieves the given property from the mesh, if any. If the 
// property is not found, this returns NULL.
void* mesh_property(mesh_t* mesh, const char* property);

// Deletes the given property from the mesh. This has no effect if the 
// property is not found.
void mesh_delete_property(mesh_t* mesh, const char* property);

// Adds a named "feature" to the mesh. A mesh either has a feature or it doesn't.
// Features can be used to make algorithmic decisions about how to perform 
// calculations on a given mesh. This function has no effect if the given 
// feature exists on the mesh.
void mesh_add_feature(mesh_t* mesh, const char* feature);

// Returns true if the mesh has the given feature, false if not.
bool mesh_has_feature(mesh_t* mesh, const char* feature);

// Deletes the given feature from the mesh if it exists. This function has 
// no effect if the mesh does not have the feature.
void mesh_delete_feature(mesh_t* mesh, const char* feature);

// Returns a newly-allocated list of indices that will define a tags for 
// cells/faces/edges/nodes with the given descriptor. If the tag already 
// exists, returns NULL.
int* mesh_create_tag(tagger_t* tagger, const char* tag, int num_indices);

// Retrieves the given tag, returning an array of indices if found (and 
// writing the number of tagged elements to num_elements), or NULL if not.
int* mesh_tag(tagger_t* tagger, const char* tag, int* num_indices);

// Returns true if the given tag exists, false if not.
bool mesh_has_tag(tagger_t* tagger, const char* tag);

// Associates a named piece of metadata (a "property") with the given tag.
// This can be used to store data related to tagged indices.
// A destructor function can be passed in to handle freeing of resources.
// If the tag is not found, this function has no effect. If the given property
// exists on the tag, it is overwritten. Returns true if the property was 
// added, false if not.
bool mesh_tag_set_property(tagger_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*));

// Retrieves the given property associated with the given tag, if any. If the 
// tag or property are not found, this returns NULL.
void* mesh_tag_property(tagger_t* tagger, const char* tag, const char* property);

// Deletes the given property from the tag. This has no effect if the tag
// or property are not found.
void mesh_tag_delete_property(tagger_t* tagger, const char* tag, const char* property);

// Renames the given tag. This has no effect if the tag is not found.
void mesh_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag);

// Deletes the given tag. This has no effect if the tag is not found.
void mesh_delete_tag(tagger_t* tagger, const char* tag);

// Allows the traversal of all mesh tags.
bool mesh_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, int* tag_size);

// Computes face areas and cell volumes for the mesh (for those that are 
// bounded).
void mesh_compute_geometry(mesh_t* mesh);

// This helper method constructs edge information from nodes, provided that 
// no edge information exists already.
void mesh_construct_edges(mesh_t* mesh);

// Returns the number of faces attached to the given cell in the mesh.
static inline int mesh_cell_num_faces(mesh_t* mesh, int cell)
{
  return mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell];
}

// Allows iteration over the faces attached to the given cell in the mesh.
// Set *pos to 0 to reset the iteration. Returns true if faces remain in 
// the cell, false otherwise.
static inline bool mesh_cell_next_face(mesh_t* mesh, int cell, int* pos, int* face)
{
  *face = mesh->cell_faces[mesh->cell_face_offsets[cell] + *pos];
  // FIXME: We don't currently use the one's complement info that 
  // FIXME: polytope provides -- we only transform indices back.
  if (*face < 0) *face = ~(*face);
  ++(*pos);
  return (*pos <= (mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell]));
}

// Allows iteration over the neighboring cells attached to the given cell in 
// the mesh, in the same order as that given by mesh_cell_next_face(). If the 
// next neighbor for a cell is non-existant, *neighbor_cell will be set to -1.
// Set *pos to 0 to reset the iteration. Returns true if the traversal over 
// all faces of the cell is not complete, false otherwise.
static inline bool mesh_cell_next_neighbor(mesh_t* mesh, int cell, int* pos, int* neighbor_cell)
{
  int face;
  bool result = mesh_cell_next_face(mesh, cell, pos, &face);
  *neighbor_cell = mesh->face_cells[2*face];
  return result;
}

// Returns the number of nodes attached to the given face in the mesh.
static inline int mesh_face_num_nodes(mesh_t* mesh, int face)
{
  return mesh->face_node_offsets[face+1] - mesh->face_node_offsets[face];
}

// Allows iteration over the nodes attached to the given face in the mesh.
// Set *pos to 0 to reset the iteration. Returns true if nodes remain in 
// the face, false otherwise.
static inline bool mesh_next_face_node(mesh_t* mesh, int face, int* pos, int* node)
{
  *node = mesh->face_nodes[mesh->face_node_offsets[face] + *pos];
  ++(*pos);
  return (*pos <= (mesh->face_node_offsets[face+1] - mesh->face_node_offsets[face]));
}

// Returns the number of edges attached to the given face in the mesh.
static inline int mesh_face_num_edges(mesh_t* mesh, int face)
{
  return mesh->face_node_offsets[face+1] - mesh->face_node_offsets[face];
}

// Allows iteration over the edges attached to the given face in the mesh.
// Set *pos to 0 to reset the iteration. Returns true if edges remain in 
// the face, false otherwise.
static inline bool mesh_next_face_edge(mesh_t* mesh, int face, int* pos, int* edge)
{
  *edge = mesh->face_edges[mesh->face_edge_offsets[face] + *pos];
  ++(*pos);
  return (*pos <= (mesh->face_edge_offsets[face+1] - mesh->face_edge_offsets[face]));
}

// Given a face within the mesh and one of its cells, returns the cell on 
// the opposite side of the face, or -1 if there is no such cell.
static inline int mesh_face_opp_cell(mesh_t* mesh, int face, int cell)
{
  return (cell == mesh->face_cells[2*face]) ? mesh->face_cells[2*face+1] 
                                            : mesh->face_cells[2*face];
}

// Returns a serializer object that can read/write meshes from/to byte arrays.
serializer_t* mesh_serializer();

// This function constructs an adjacency graph expressing the connectivity of 
// the cells of the given mesh.
adj_graph_t* graph_from_mesh_cells(mesh_t* mesh);

#endif

