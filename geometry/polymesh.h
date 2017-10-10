// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYMESH_H
#define POLYMEC_POLYMESH_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/tagger.h"
#include "core/exchanger.h"
#include "core/adj_graph.h"
#include "core/serializer.h"

// Mesh centerings.
typedef enum
{
  POLYMESH_NODE,
  POLYMESH_EDGE,
  POLYMESH_FACE,
  POLYMESH_CELL
} polymesh_centering_t;

// Mesh features.
extern const char* POLYMESH_IS_TETRAHEDRAL; // indicates that a polymesh is tetrahedral.

typedef struct polymesh_storage_t polymesh_storage_t;

// This data type represents an unstructured polyhedral mesh, consisting of 
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
  polymesh_storage_t* storage;
} polymesh_t;

// Construct a new polymesh with the given number of cells, ghost cells, 
// faces, and nodes. This function does not provide any description
// of the mesh's topology and is only useful in the construction of mesh 
// generation algorithms. 
// NOTE: edges are constructed with polymesh_construct_edges().
polymesh_t* polymesh_new(MPI_Comm comm, int num_cells, int num_ghost_cells, 
                         int num_faces, int num_nodes);

// Construct a new polymesh with a single type of polytope.
// This function does not set up connectivity, but initializes its metadata 
// according to the prescribed number of faces per cell and nodes per face.
// No edge connectivity is set up.
polymesh_t* polymesh_new_with_cell_type(MPI_Comm comm, int num_cells, 
                                        int num_ghost_cells, int num_faces, 
                                        int num_nodes, int num_faces_per_cell,
                                        int num_nodes_per_face);

// Destroys the given polymesh.
void polymesh_free(polymesh_t* mesh);

// Verifies the topological correctness of the polymesh, calling the given 
// (variadic) handler function with a formatted string containing a 
// description of any problems encountered. If the topology of the mesh is 
// correct, this function returns true and the handler function is not called.
// Otherwise, the function returns false.
bool polymesh_verify_topology(polymesh_t* mesh, 
                              void (*handler)(const char* format, ...));

// Returns an exact copy of the given polymesh.
polymesh_t* polymesh_clone(polymesh_t* mesh);

// Associates a named piece of metadata (a "property") with the polymesh itself.
// This can be used to store information about (for example) how the mesh 
// was generated, which can sometimes be useful. A serializer can 
// be given so that any partitioning or repartitioning of the mesh can 
// preserve this property on subdomains. If the given property exists on the 
// mesh, it is overwritten.
void polymesh_set_property(polymesh_t* mesh, 
                           const char* property, 
                           void* data, 
                           serializer_t* serializer);

// Retrieves the given property from the polymesh, if any. If the 
// property is not found, this returns NULL.
void* polymesh_property(polymesh_t* mesh, const char* property);

// Deletes the given property from the polymesh. This has no effect if the 
// property is not found.
void polymesh_delete_property(polymesh_t* mesh, const char* property);

// Allows the traversal of polymesh properties. Set *pos to 0 to reset the 
// iteration.
bool polymesh_next_property(polymesh_t* mesh, int* pos, 
                            char** prop_name, void** prop_data, 
                            serializer_t** prop_serializer);

// Returns an exchanger object that can be used to perform parallel exchanges
// on cell-centered polymesh data. In serial configurations, this exchanger holds 
// no data and exchanges have no effect.
exchanger_t* polymesh_exchanger(polymesh_t* mesh);

// Sets the exchanger to be used by the polymesh, replacing any existing exchanger.
void polymesh_set_exchanger(polymesh_t* mesh, exchanger_t* ex);

// Adds a named "feature" to the polymesh. A mesh either has a feature or it doesn't.
// Features can be used to make algorithmic decisions about how to perform 
// calculations on a given mesh. This function has no effect if the given 
// feature exists on the mesh.
void polymesh_add_feature(polymesh_t* mesh, const char* feature);

// Returns true if the polymesh has the given feature, false if not.
bool polymesh_has_feature(polymesh_t* mesh, const char* feature);

// Deletes the given feature from the polymesh if it exists. This function has 
// no effect if the mesh does not have the feature.
void polymesh_delete_feature(polymesh_t* mesh, const char* feature);

// Returns a newly-allocated list of indices that will define a tags for 
// cells/faces/edges/nodes with the given descriptor. If the tag already 
// exists, returns NULL.
int* polymesh_create_tag(tagger_t* tagger, const char* tag, size_t num_indices);

// Retrieves the given tag, returning an array of indices if found (and 
// writing the number of tagged elements to num_elements), or NULL if not.
int* polymesh_tag(tagger_t* tagger, const char* tag, size_t* num_indices);

// Returns true if the given tag exists, false if not.
bool polymesh_has_tag(tagger_t* tagger, const char* tag);

// Associates a named piece of metadata (a "property") with the given tag.
// This can be used to store data related to tagged indices.
// A serializer should be provided for the property's data.
// If the tag is not found, this function has no effect. If the given property
// exists on the tag, it is overwritten. Returns true if the property was 
// added, false if not.
bool polymesh_tag_set_property(tagger_t* tagger, const char* tag, const char* property, void* data, serializer_t* serializer);

// Retrieves the given property associated with the given tag, if any. If the 
// tag or property are not found, this returns NULL.
void* polymesh_tag_property(tagger_t* tagger, const char* tag, const char* property);

// Deletes the given property from the tag. This has no effect if the tag
// or property are not found.
void polymesh_tag_delete_property(tagger_t* tagger, const char* tag, const char* property);

// Allows the traversal of properties on a tag. Set *pos to 0 to reset the 
// iteration.
bool polymesh_tag_next_property(tagger_t* tagger, const char* tag, int* pos, 
                                char** prop_name, void** prop_data, 
                                serializer_t** prop_serializer);

// Renames the given tag. This has no effect if the tag is not found.
void polymesh_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag);

// Deletes the given tag. This has no effect if the tag is not found.
void polymesh_delete_tag(tagger_t* tagger, const char* tag);

// Allows the traversal of all polymesh tags.
bool polymesh_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, size_t* tag_size);

// Computes face areas and cell volumes for the polymesh (for those that are 
// bounded).
void polymesh_compute_geometry(polymesh_t* mesh);

// This helper method makes sure that sufficient storage is reserved for 
// cell-face and face->node connectivity. This must be called after 
// the mesh->cell_face_offsets and mesh->face_node_offsets arrays have been 
// properly initialized. It does *NOT* allocate edge information -- use 
// polymesh_construct_edges() to build face->edge and edge->node connectivity.
void polymesh_reserve_connectivity_storage(polymesh_t* mesh);

// This helper method constructs face->edge and edge->node connectivity, 
// provided that no edge information exists already.
void polymesh_construct_edges(polymesh_t* mesh);

// Returns the number of faces attached to the given cell in the mesh.
static inline int polymesh_cell_num_faces(polymesh_t* mesh, int cell)
{
  return mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell];
}

// Allows iteration over the faces attached to the given cell in the polymesh.
// Set *pos to 0 to reset the iteration. Returns true if faces remain in 
// the cell, false otherwise. NOTE: the local index of the face within the 
// cell is *pos - 1 after the call. This method always returns a non-negative 
// face index.
static inline bool polymesh_cell_next_face(polymesh_t* mesh, int cell, int* pos, int* face)
{
  *face = mesh->cell_faces[mesh->cell_face_offsets[cell] + *pos];
  if (*face < 0) *face = ~(*face);
  ++(*pos);
  return (*pos <= (mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell]));
}

// Allows iteration over the oriented faces attached to the given cell in the 
// polymesh. Set *pos to 0 to reset the iteration. Returns true if faces remain in 
// the cell, false otherwise. NOTE: the local index of the face within the 
// cell is *pos - 1 after the call. This method returns a non-negative face index 
// if the nodes in the face are to be traversed in order, or the (negative)
// two's complement to the actual face index if its nodes are to be traversed 
// in reverse order.
static inline bool polymesh_cell_next_oriented_face(polymesh_t* mesh, int cell, int* pos, int* face)
{
  *face = mesh->cell_faces[mesh->cell_face_offsets[cell] + *pos];
  ++(*pos);
  return (*pos <= (mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell]));
}

// Allows iteration over the neighboring cells attached to the given cell in 
// the polymesh, in the same order as that given by polymesh_cell_next_face(). If the 
// next neighbor for a cell is non-existent, *neighbor_cell will be set to -1.
// Set *pos to 0 to reset the iteration. Returns true if the traversal over 
// all faces of the cell is not complete, false otherwise. NOTE: the local 
// index of the face separating the cells is *pos - 1 after the call.
static inline bool polymesh_cell_next_neighbor(polymesh_t* mesh, int cell, int* pos, int* neighbor_cell)
{
  int face;
  bool result = polymesh_cell_next_face(mesh, cell, pos, &face);
  if (mesh->face_cells[2*face] == cell)
    *neighbor_cell = mesh->face_cells[2*face+1];
  else
    *neighbor_cell = mesh->face_cells[2*face];
  return result;
}

// This returns the index of the face shared by cell and neighbor_cell if the two 
// cells share a face, -1 otherwise.
static inline int polymesh_cell_face_for_neighbor(polymesh_t* mesh, int cell, int neighbor_cell)
{
  int pos = 0, face;
  while (polymesh_cell_next_face(mesh, cell, &pos, &face))
  {
    if ((mesh->face_cells[2*face] == neighbor_cell) || (mesh->face_cells[2*face+1] == neighbor_cell))
      return face;
  }
  return -1;
}

// Returns true if cell1 and cell2 are neighbors that share a face, 
// false otherwise.
static inline bool polymesh_cells_are_neighbors(polymesh_t* mesh, int cell1, int cell2)
{
  return (polymesh_cell_face_for_neighbor(mesh, cell1, cell2) != -1);
}

// Returns the number of nodes attached to the given face in the mesh.
static inline int polymesh_face_num_nodes(polymesh_t* mesh, int face)
{
  return mesh->face_node_offsets[face+1] - mesh->face_node_offsets[face];
}

// Allows iteration over the nodes attached to the given face in the polymesh.
// Set *pos to 0 to reset the iteration. Returns true if nodes remain in 
// the face, false otherwise. NOTE: the local index of the node within the 
// face is *pos - 1 after the call. If face is the (negative) two's complement 
// of the actual face index, the nodes of the face will be traversed in reverse
// order.
static inline bool polymesh_face_next_node(polymesh_t* mesh, int face, int* pos, int* node)
{
  int actual_face;
  if (face >= 0)
  {
    actual_face = face;
    *node = mesh->face_nodes[mesh->face_node_offsets[actual_face] + *pos];
  }
  else
  {
    actual_face = ~face;
    // We have to take care not to step off the beginning of the face_nodes array.
    int offset = mesh->face_node_offsets[actual_face+1] - *pos - 1;
    if (offset >= 0)
      *node = mesh->face_nodes[mesh->face_node_offsets[actual_face+1] - *pos - 1];
  }
  ++(*pos);
  return (*pos <= (mesh->face_node_offsets[actual_face+1] - mesh->face_node_offsets[actual_face]));
}

// Returns the number of edges attached to the given face in the mesh.
static inline int polymesh_face_num_edges(polymesh_t* mesh, int face)
{
  return mesh->face_node_offsets[face+1] - mesh->face_node_offsets[face];
}

// Allows iteration over the edges attached to the given face in the mesh.
// Set *pos to 0 to reset the iteration. Returns true if edges remain in 
// the face, false otherwise. NOTE: the local index of the edge within the 
// face is *pos - 1 after the call. If face is the (negative) two's complement 
// of the actual face index, the edges of the face will be traversed in reverse
// order.
static inline bool polymesh_face_next_edge(polymesh_t* mesh, int face, int* pos, int* edge)
{
  int actual_face;
  if (face >= 0)
  {
    actual_face = face;
    *edge = mesh->face_edges[mesh->face_edge_offsets[actual_face] + *pos];
  }
  else
  {
    actual_face = ~face;
    // We have to take care not to step off the beginning of the face_edges array.
    int offset = mesh->face_edge_offsets[actual_face+1] - *pos - 1;
    if (offset >= 0)
      *edge = mesh->face_edges[mesh->face_edge_offsets[actual_face+1] - *pos - 1];
  }
  ++(*pos);
  return (*pos <= (mesh->face_edge_offsets[actual_face+1] - mesh->face_edge_offsets[actual_face]));
}

// This returns the index of the edge shared by face and neighbor_face if the two 
// faces share a edge, -1 otherwise. A non-negative face index must be given.
// This function can be somewhat costly, since it requires a linear search through 
// the edges of the two faces.
static inline int polymesh_face_edge_for_neighbor(polymesh_t* mesh, int face, int neighbor_face)
{
  for (int e1 = mesh->face_edge_offsets[face]; e1 < mesh->face_edge_offsets[face+1]; ++e1)
  {
    int edge1 = mesh->face_edges[e1];
    for (int e2 = mesh->face_edge_offsets[neighbor_face]; e2 < mesh->face_edge_offsets[neighbor_face+1]; ++e2)
    {
      int edge2 = mesh->face_edges[e2];
      if (edge2 == edge1)
        return edge2;
    }
  }
  return -1;
}

// Returns the "first" cell attached to a face. A well-formed face has at least one 
// cell with a non-negative index.
static inline int polymesh_face_cell1(polymesh_t* mesh, int face)
{
  return mesh->face_cells[2*face];
}

// Returns the "second" cell attached to a face. If the face is only attached to 
// one cell, the second cell is -1.
static inline int polymesh_face_cell2(polymesh_t* mesh, int face)
{
  return mesh->face_cells[2*face+1];
}

// Given a face within the polymesh and one of its cells, returns the cell on 
// the opposite side of the face, or -1 if there is no such cell.
static inline int polymesh_face_opp_cell(polymesh_t* mesh, int face, int cell)
{
  return (cell == mesh->face_cells[2*face]) ? mesh->face_cells[2*face+1] 
                                            : mesh->face_cells[2*face];
}

// Returns true if the given face in the polymesh abuts an external boundary, 
// false if it has an opposite cell.
static inline bool polymesh_face_is_external(polymesh_t* mesh, int face)
{
  return (mesh->face_cells[2*face+1] == -1);
}

// Returns a serializer object that can read/write polymeshes from/to byte arrays.
serializer_t* polymesh_serializer(void);

// 1-value face exchanger: creates and returns a newly-allocated exchanger that allows the 
// exchange of a UNIQUE face-related value from local to remote processes. The array to be 
// exchanged should be of length stride * mesh->num_faces, and the data is 
// laid out thus:
// data[face*stride + s] is the sth value of the data associated with 
// the local face. The process with the lowest rank on which the face appears
// owns that face. Communication is required to construct this face exchanger.
exchanger_t* polymesh_1v_face_exchanger_new(polymesh_t* mesh);

// 2-value face exchanger: creates and returns a newly-allocated exchanger that allows the 
// exchange of BOTH face-related values from local to remote processes. The array to be 
// exchanged should be of length 2 * stride * mesh->num_faces, and the data is 
// laid out thus:
// data[(2*face+c)*stride + s] is the sth value of the data associated with 
// the cth (0th or 1st) value of the (local) face. Local values are associated 
// with cell 0 of local faces on parallel domain boundaries, and remote values
// are associated with cell 1.
// No communication is required to construct this face exchanger.
exchanger_t* polymesh_2v_face_exchanger_new(polymesh_t* mesh);

// 1-value node exchanger: creates and returns a newly-allocated exchanger that allows the 
// population of node-centered arrays with a unique value for each node shared between 
// processes. The process on which a node is represented with the lowest rank owns that node.
// NOTE: Currently, in order to construct a node exchanger, one must begin with a mesh whose node positions 
// are geometrically consistent in the sense that each node on a domain is closest to or 
// colocated with exactly one node on one or more other domains. If this requirement is not 
// met, the node exchanger cannot be reliably created.
// Communication is required to construct such a node exchanger.
exchanger_t* polymesh_1v_node_exchanger_new(polymesh_t* mesh);

// n-value node exchanger: creates and returns a newly-allocated exchanger that allows the 
// exchange of multiple node-related values from local to remote processes, and populates 
// the given node_offsets array (of length mesh->num_nodes + 1) with offsets that describe 
// the layout of arrays to be arranged. The array to be exchanged should be of length 
// node_offsets[mesh->num_nodes] * stride, and is laid out as follows for an array A:
// - the local nodal value for node n is located at A[node_offsets[n-1]].
// - the values of A at indices node_offsets[n-1]+1 through node_offsets[n] - 1 contain the 
//   remote values of the node n as they exist on other processes.
// NOTE: Currently, in order to construct a node exchanger, one must begin with a mesh whose node positions 
// are geometrically consistent in the sense that each node on a domain is closest to or 
// colocated with exactly one node on one or more other domains. If this requirement is not 
// met, the node exchanger cannot be reliably created.
// Communication is required to construct such a node exchanger.
exchanger_t* polymesh_nv_node_exchanger_new(polymesh_t* mesh, int* node_offsets);

// This function constructs an adjacency graph expressing the connectivity of 
// the cells of the given polymesh.
adj_graph_t* graph_from_polymesh_cells(polymesh_t* mesh);

#endif

