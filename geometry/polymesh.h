// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYMESH_H
#define POLYMEC_POLYMESH_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/adj_graph.h"
#include "core/exchanger.h"
#include "core/serializer.h"
#include "geometry/tagger.h"

/// \addtogroup geometry geometry
///@{

/// \enum polymesh_centering_t
/// Centerings for polyhedral meshes (polymeshes).
typedef enum
{
  POLYMESH_NODE,
  POLYMESH_EDGE,
  POLYMESH_FACE,
  POLYMESH_CELL
} polymesh_centering_t;

typedef struct polymesh_storage_t polymesh_storage_t;

/// \class polymesh
/// This type represents an unstructured polyhedral mesh consisting of
/// stationary cells and the faces connecting them.
struct polymesh_t
{
  /// MPI communicator.
  MPI_Comm comm;

  /// The total number of (locally-owned) cells in the mesh.
  int num_cells;
  /// The number of ghost cells on the local domain in the mesh.
  int num_ghost_cells;
  /// The offsets of the sets of faces attached to cells, stored in CRS format.
  int* cell_face_offsets;
  /// The indices of faces attached to cells, stored in CRS format.
  int* cell_faces;

  /// The total number of local faces in the mesh.
  int num_faces;
  /// The offsets of the sets of nodes attached to faces, stored in CRS format.
  int* face_node_offsets;
  /// The indices of nodes attached to faces, stored in CRS format.
  int* face_nodes;

  /// The offsets of the sets of edges attached to faces, stored in CRS format.
  int* face_edge_offsets;
  /// The indices of edges attached to faces, stored in CRS format.
  int* face_edges;

  /// The cells attached to the faces in the mesh. Each face has 2 cells,
  /// so the first cell of the ith face is face_cells[2*i] and the second
  /// is face_cells[2*i+1].
  int* face_cells;

  /// The total number of local edges in the mesh.
  int num_edges;

  /// The nodes of the edges in the mesh. Each edge has 2 nodes, so the
  /// the first node of the ith edge is edge_nodes[2*i] and the second is
  /// edge_nodes[2*i+1].
  int* edge_nodes;

  /// The total number of local nodes in the mesh.
  int num_nodes;
  /// Coordinates of the mesh nodes, indexed from 0 to N-1.
  point_t* nodes;

  /// Geometry information.
  real_t* cell_volumes;
  point_t* cell_centers;
  point_t* face_centers;
  real_t* face_areas;
  vector_t* face_normals;

  ///@{
  /// Mesh tagging mechanisms.
  tagger_t* cell_tags;
  tagger_t* face_tags;
  tagger_t* edge_tags;
  tagger_t* node_tags;
  ///@}

  /// Mesh storage information -- used internally.
  polymesh_storage_t* storage;
};
typedef struct polymesh_t polymesh_t;

/// Constructs a new polymesh with the given number of cells, ghost cells,
/// faces, and nodes. This function does not provide any description
/// of the mesh's topology and is only useful in the construction of mesh
/// generation algorithms.
/// \note edges are constructed with polymesh_construct_edges().
/// \memberof polymesh
polymesh_t* polymesh_new(MPI_Comm comm, int num_cells, int num_ghost_cells,
                         int num_faces, int num_nodes);

/// Constructs a new polymesh with a single type of polytope.
/// This function does not set up connectivity, but initializes its metadata
/// according to the prescribed number of faces per cell and nodes per face.
/// No edge connectivity is set up.
/// \memberof polymesh
polymesh_t* polymesh_new_with_cell_type(MPI_Comm comm, int num_cells,
                                        int num_ghost_cells, int num_faces,
                                        int num_nodes, int num_faces_per_cell,
                                        int num_nodes_per_face);

/// Destroys the given polymesh.
/// \memberof polymesh
void polymesh_free(polymesh_t* mesh);

/// Returns true if the given mesh is topologically correct, false if not.
/// \param [out] reason If the mesh is not topologically correct and this is
///                     non-NULL, it stores an internal pointer to a string
///                     that explains why.
/// \memberof polymesh
bool polymesh_is_valid(polymesh_t* mesh, char** reason);

/// Returns an exact copy of the given polymesh.
/// \memberof polymesh
polymesh_t* polymesh_clone(polymesh_t* mesh);

/// Returns an exchanger object that can be used to perform parallel exchanges
/// on polymesh fields with the given centering. In serial configurations,
/// this exchanger holds no data, and exchanges have no effect.
/// \param [in] centering The centering of the data handled by this exchanger.
/// \memberof polymesh
exchanger_t* polymesh_exchanger(polymesh_t* mesh,
                                polymesh_centering_t centering);

/// Returns a newly-allocated list of indices that will define a tags for
/// cells/faces/edges/nodes with the given descriptor. If the tag already
/// exists, returns NULL.
/// \memberof polymesh
int* polymesh_create_tag(tagger_t* tagger, const char* tag, size_t num_indices);

/// Retrieves the given tag, returning an array of indices if found (and
/// writing the number of tagged elements to num_elements), or NULL if not.
/// \memberof polymesh
int* polymesh_tag(tagger_t* tagger, const char* tag, size_t* num_indices);

/// Returns true if the given tag exists, false if not.
/// \memberof polymesh
bool polymesh_has_tag(tagger_t* tagger, const char* tag);

/// Renames the given tag. This has no effect if the tag is not found.
/// \memberof polymesh
void polymesh_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag);

/// Deletes the given tag. This has no effect if the tag is not found.
/// \memberof polymesh
void polymesh_delete_tag(tagger_t* tagger, const char* tag);

/// Allows the traversal of all polymesh tags.
/// \memberof polymesh
bool polymesh_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, size_t* tag_size);

/// Computes face areas and cell volumes for the polymesh (for those that are
/// bounded).
/// \memberof polymesh
void polymesh_compute_geometry(polymesh_t* mesh);

/// This helper method makes sure that sufficient storage is reserved for
/// cell-face and face->node connectivity. This must be called after
/// the mesh->cell_face_offsets and mesh->face_node_offsets arrays have been
/// properly initialized. It does *NOT* allocate edge information -- use
/// polymesh_construct_edges() to build face->edge and edge->node connectivity.
/// \memberof polymesh
void polymesh_reserve_connectivity_storage(polymesh_t* mesh);

/// This helper method constructs face->edge and edge->node connectivity,
/// provided that no edge information exists already.
/// \memberof polymesh
void polymesh_construct_edges(polymesh_t* mesh);

/// Returns the number of faces attached to the given cell in the mesh.
/// \memberof polymesh
static inline int polymesh_cell_num_faces(polymesh_t* mesh, int cell)
{
  return mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell];
}

/// Retrieves the faces attached to the given cell in the polymesh.
/// \param [out] faces An array large enough to store all the faces for the given cell.
/// \memberof polymesh
static inline void polymesh_cell_get_faces(polymesh_t* mesh, int cell, int* faces)
{
  int start = mesh->cell_face_offsets[cell];
  int end = mesh->cell_face_offsets[cell+1];
  for (int f = start; f < end; ++f)
    faces[f-start] = (mesh->cell_faces[f] < 0) ? ~(mesh->cell_faces[f]) : mesh->cell_faces[f];
}

/// Allows iteration over the faces attached to the given cell in the polymesh.
/// Set *pos to 0 to reset the iteration. Returns true if faces remain in
/// the cell, false otherwise. NOTE: the local index of the face within the
/// cell is *pos - 1 after the call. This method always returns a non-negative
/// face index.
/// \memberof polymesh
static inline bool polymesh_cell_next_face(polymesh_t* mesh, int cell, int* pos, int* face)
{
  *face = mesh->cell_faces[mesh->cell_face_offsets[cell] + *pos];
  if (*face < 0) *face = ~(*face);
  ++(*pos);
  return (*pos <= (mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell]));
}

/// Retrieves the oriented faces attached to the given cell in the polymesh.
/// \param [out] faces An array large enough to store all the faces for the given
///                    cell. If the face is negative, it points inward relative
///                    to the cell, and its 1's complement is its index.
/// \memberof polymesh
static inline void polymesh_cell_get_oriented_faces(polymesh_t* mesh, int cell, int* faces)
{
  int start = mesh->cell_face_offsets[cell];
  int end = mesh->cell_face_offsets[cell+1];
  for (int f = start; f < end; ++f)
    faces[f-start] = mesh->cell_faces[f];
}

/// Allows iteration over the oriented faces attached to the given cell in the
/// polymesh.
/// \param [in] cell The cell whose faces are traversed.
/// \param [in,out] pos Controls the iteration. Set *pos to 0 to reset.
/// \param [out] face Stores the next face in the traversal. Stores a non-negative face index
///                   if the nodes in the face are to be traversed in order, or the (negative)
///                   one's complement to the actual face index if its nodes are to be traversed
///                   in reverse order.
/// \returns true if faces remain in the cell, false otherwise.
/// \note The local index of the face within the cell is *pos - 1 after the call.
/// \memberof polymesh
static inline bool polymesh_cell_next_oriented_face(polymesh_t* mesh, int cell, int* pos, int* face)
{
  *face = mesh->cell_faces[mesh->cell_face_offsets[cell] + *pos];
  ++(*pos);
  return (*pos <= (mesh->cell_face_offsets[cell+1] - mesh->cell_face_offsets[cell]));
}

/// Allows iteration over the neighboring cells attached to the given cell in
/// the polymesh, in the same order as that given by polymesh_cell_next_face(). If the
/// next neighbor for a cell is non-existent, *neighbor_cell will be set to -1.
/// Set *pos to 0 to reset the iteration. Returns true if the traversal over
/// all faces of the cell is not complete, false otherwise. NOTE: the local
/// index of the face separating the cells is *pos - 1 after the call.
/// \memberof polymesh
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

/// This returns the index of the face shared by cell and neighbor_cell if the two
/// cells share a face, -1 otherwise.
/// \memberof polymesh
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

/// Returns true if cell1 and cell2 are neighbors that share a face,
/// false otherwise.
/// \memberof polymesh
static inline bool polymesh_cells_are_neighbors(polymesh_t* mesh, int cell1, int cell2)
{
  return (polymesh_cell_face_for_neighbor(mesh, cell1, cell2) != -1);
}

/// Returns the number of nodes attached to the given face in the mesh.
/// \memberof polymesh
static inline int polymesh_face_num_nodes(polymesh_t* mesh, int face)
{
  return mesh->face_node_offsets[face+1] - mesh->face_node_offsets[face];
}

/// Retrieves the nodes attached to the given face in the polymesh.
/// \param [out] nodes An array large enough to store all the nodes for the given
///                    face.
/// \memberof polymesh
static inline void polymesh_face_get_nodes(polymesh_t* mesh, int face, int* nodes)
{
  if (face < 0)
    face = ~face;
  int start = mesh->face_node_offsets[face];
  int end = mesh->face_node_offsets[face+1];
  for (int n = start; n < end; ++n)
    nodes[n-start] = mesh->face_nodes[n];
}

/// Allows iteration over the nodes attached to the given face in the polymesh.
/// Set *pos to 0 to reset the iteration. Returns true if nodes remain in
/// the face, false otherwise. NOTE: the local index of the node within the
/// face is *pos - 1 after the call. If face is the (negative) two's complement
/// of the actual face index, the nodes of the face will be traversed in reverse
/// order.
/// \memberof polymesh
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

/// Returns the number of edges attached to the given face in the mesh.
/// \memberof polymesh
static inline int polymesh_face_num_edges(polymesh_t* mesh, int face)
{
  return mesh->face_edge_offsets[face+1] - mesh->face_edge_offsets[face];
}

/// Retrieves the edges attached to the given face in the polymesh.
/// \param [out] edges An array large enough to store all the edges for the given
///                    face.
/// \memberof polymesh
static inline void polymesh_face_get_edges(polymesh_t* mesh, int face, int* edges)
{
  if (face < 0)
    face = ~face;
  int start = mesh->face_edge_offsets[face];
  int end = mesh->face_edge_offsets[face+1];
  for (int e = start; e < end; ++e)
    edges[e-start] = mesh->face_edges[e];
}

/// Allows iteration over the edges attached to the given face in the mesh.
/// Set *pos to 0 to reset the iteration. Returns true if edges remain in
/// the face, false otherwise. NOTE: the local index of the edge within the
/// face is *pos - 1 after the call. If face is the (negative) two's complement
/// of the actual face index, the edges of the face will be traversed in reverse
/// order.
/// \memberof polymesh
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

/// This returns the index of the edge shared by face and neighbor_face if the two
/// faces share a edge, -1 otherwise. A non-negative face index must be given.
/// This function can be somewhat costly, since it requires a linear search through
/// the edges of the two faces.
/// \memberof polymesh
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

/// Returns the "first" cell attached to a face. A well-formed face has at least one
/// cell with a non-negative index.
/// \memberof polymesh
static inline int polymesh_face_cell1(polymesh_t* mesh, int face)
{
  return mesh->face_cells[2*face];
}

/// Returns the "second" cell attached to a face. If the face is only attached to
/// one cell, the second cell is -1.
/// \memberof polymesh
static inline int polymesh_face_cell2(polymesh_t* mesh, int face)
{
  return mesh->face_cells[2*face+1];
}

/// Given a face within the polymesh and one of its cells, returns the cell on
/// the opposite side of the face, or -1 if there is no such cell.
/// \memberof polymesh
static inline int polymesh_face_opp_cell(polymesh_t* mesh, int face, int cell)
{
  return (cell == mesh->face_cells[2*face]) ? mesh->face_cells[2*face+1]
                                            : mesh->face_cells[2*face];
}

/// Returns true if the given face in the polymesh abuts an external boundary,
/// false if it has an opposite cell.
/// \memberof polymesh
static inline bool polymesh_face_is_external(polymesh_t* mesh, int face)
{
  return (mesh->face_cells[2*face+1] == -1);
}

/// Returns a serializer object that can read/write polymeshes from/to byte arrays.
/// \memberof polymesh
serializer_t* polymesh_serializer(void);

/// This function constructs an adjacency graph expressing the connectivity of
/// the cells of the given polymesh.
/// \relates polymesh
adj_graph_t* graph_from_polymesh_cells(polymesh_t* mesh);

///@}

#endif

