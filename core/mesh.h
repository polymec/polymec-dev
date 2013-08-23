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

#ifndef POLYMEC_MESH_H
#define POLYMEC_MESH_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/sp_func.h"
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
typedef struct mesh_tags_t mesh_tags_t;

typedef struct face_t face_t;

// A node in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} node_t;

// This helper function computes the displacement vector pointing from 
// node1 to node2, storing it in displacement.
static inline void node_displacement(node_t* node1, node_t* node2, vector_t* displacement)
{
  displacement->x = node2->x - node1->x;
  displacement->y = node2->y - node1->y;
  displacement->z = node2->z - node1->z;
}

// An edge, which is bounded by two nodes.
typedef struct 
{
  // Nodes that bound this edge. node1 is always the node with the lower index.
  node_t* node1;
  node_t* node2;
  // Length of the edge.
  double length;
} edge_t;

// A full-featured cell. This cell is bounded by a set of faces.
typedef struct
{
  // Faces bounding this cell. If this cell is a ghost cell, the faces 
  // do not bound the entire cell--only those that connect this cell 
  // with interior cells are included.
  face_t** faces;
  // Number of faces attached to this cell.
  int num_faces;
  // Cell center position.
  point_t center;
  // Cell volume.
  double volume;
} cell_t;

// Returns the edge within the given cell whose nodes are node1 and node2,
// or NULL if the cell contains no such edge. The order of the nodes in 
// the edge does not matter.
edge_t* cell_find_edge_with_nodes(cell_t* cell, node_t* node1, node_t* node2);

// A full-featured face. This face connects exactly two cells and is 
// bounded by edges.
struct face_t
{
  // The cells connected by this face. cell1 is always the face with 
  // the lower index in the mesh.
  cell_t* cell1;
  cell_t* cell2;
  // Edges that bound this face.
  edge_t** edges;
  // The number of edges bounding this face.
  int num_edges;
  // Face center.
  point_t center;
  // Face normal vector (scaled by area).
  vector_t normal;
};

// This function returns the cell opposite the given cell through the 
// given face.
static inline cell_t* face_opp_cell(face_t* face, cell_t* cell)
{
  return (face->cell1 == cell) ? face->cell2 : face->cell1;
}

// This function fills the given array with the nodes of its constituent 
// edges, ordered in one of the two possible traversal orders. The number of 
// nodes is less than or equal to the number of edges. This fills the array 
// of nodes, assuming that it is properly sized for the face.
void face_get_nodes(face_t* face, node_t** nodes);
 
typedef struct mesh_storage_t mesh_storage_t;

// This data type represents an unstructured mesh, consisting of 
// stationary cells and the faces connecting them.
typedef struct 
{
  // Mesh cells, indexed from 0 to C-1.
  cell_t* cells;
  // Total number of (locally-owned) cells in the mesh.
  int num_cells;
  // Total number of ghost cells in the mesh.
  int num_ghost_cells;
  // Mesh faces, indexed from 0 to F-1.
  face_t* faces;
  // Total number of faces in the mesh.
  int num_faces;
  // Mesh edges, indexed from 0 to E-1.
  edge_t* edges;
  // Total number of edges in the mesh.
  int num_edges;
  // Mesh nodes, indexed from 0 to N-1.
  node_t* nodes;
  // Total number of nodes in the mesh.
  int num_nodes;

  // Mesh tagging mechanisms.
  mesh_tags_t* cell_tags;
  mesh_tags_t* face_tags;
  mesh_tags_t* edge_tags;
  mesh_tags_t* node_tags;

  // Mesh storage information -- used internally.
  ARENA* arena;
  bool close_arena;
  mesh_storage_t* storage;
} mesh_t;

// Construct a new mesh with the given number of cells, ghost cells, 
// faces, edges, and nodes. This function does not provide any description
// of the mesh's topology and is only useful in the construction of mesh 
// generation algorithms.
mesh_t* mesh_new(int num_cells, int num_ghost_cells, int num_faces,
                 int num_edges, int num_nodes);

// Construct a new mesh, using the given arena for memory allocations.
mesh_t* mesh_new_with_arena(ARENA* arena, int num_cells, int num_ghost_cells, int num_faces,
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

// Returns a newly-allocated list of indices that will define a tags for 
// cells/faces/edges/nodes with the given descriptor. If the tag already 
// exists, returns NULL.
int* mesh_create_tag(mesh_tags_t* tagger, const char* tag, int num_indices);

// Retrieves the given tag, returning an array of indices if found (and 
// writing the number of tagged elements to num_elements), or NULL if not.
int* mesh_tag(mesh_tags_t* tagger, const char* tag, int* num_indices);

// Returns true if the given tag exists, false if not.
bool mesh_has_tag(mesh_tags_t* tagger, const char* tag);

// Associates a named piece of metadata (a "property") with the given tag.
// This can be used to store data related to tagged indices.
// A destructor function can be passed in to handle freeing of resources.
// If the tag is not found, this function has no effect. If the given property
// exists on the tag, it is overwritten. Returns true if the property was 
// added, false if not.
bool mesh_tag_set_property(mesh_tags_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*));

// Retrieves the given property associated with the given tag, if any. If the 
// tag or property are not found, this returns NULL.
void* mesh_tag_property(mesh_tags_t* tagger, const char* tag, const char* property);

// Deletes the given property from the tag. This has no effect if the tag
// or property are not found.
void mesh_tag_delete_property(mesh_tags_t* tagger, const char* tag, const char* property);

// Renames the given tag. This has no effect if the tag is not found.
void mesh_rename_tag(mesh_tags_t* tagger, const char* old_tag, const char* new_tag);

// Deletes the given tag. This has no effect if the tag is not found.
void mesh_delete_tag(mesh_tags_t* tagger, const char* tag);

// Maps the mesh from its existing coordinates to a set of coordinates
// defined by the (vector-valued) tranformation function.
void mesh_map(mesh_t* mesh, sp_func_t* mapping);

// Computes face areas and cell volumes for the mesh (for those that are 
// bounded).
void mesh_compute_geometry(mesh_t* mesh);

// Performs some basic mesh consistency checks and calls polymec_error 
// if any inconsistencies are found. NOTE that you shouldn't call this 
// on a mesh with unbounded cells, since those cells are not "consistent" 
// in the usual sense.
void mesh_check_consistency(mesh_t* mesh);

// Computes the cell volume and center for a bounded cell.
void cell_compute_geometry(cell_t* cell);

// This prints a text representation of the given node to the given FILE.
// It needs to mesh to convey its topology.
void node_fprintf(node_t* node, mesh_t* mesh, FILE* stream);

// This prints a text representation of the given edge to the given FILE.
// It needs to mesh to convey its topology.
void edge_fprintf(edge_t* edge, mesh_t* mesh, FILE* stream);

// This prints a text representation of the given face to the given FILE.
// It needs to mesh to convey its topology.
void face_fprintf(face_t* face, mesh_t* mesh, FILE* stream);

// This prints a text representation of the given cell to the given FILE.
// It needs to mesh to convey its topology.
void cell_fprintf(cell_t* cell, mesh_t* mesh, FILE* stream);

#endif

