#ifndef ARBI_MESH_H
#define ARBI_MESH_H

#include "core/arbi.h"
#include "core/point.h"
#include "arena/proto.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  // Face area.
  double area;
  // Face center.
  point_t center;
};

static inline cell_t* face_opp_cell(face_t* face, cell_t* cell)
{
  return (face->cell1 == cell) ? face->cell2 : face->cell1;
}
 
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

// Deletes the given tag. This has no effect if the tag is not found.
void mesh_delete_tag(mesh_tags_t* tagger, const char* tag);

#ifdef __cplusplus
}
#endif

#endif

