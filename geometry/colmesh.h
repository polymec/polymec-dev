// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_COLMESH_H
#define POLYMEC_COLMESH_H

#include "core/point.h"
#include "core/exchanger.h"
#include "geometry/planar_polymesh.h"
#include "geometry/polygon.h"

/// \addtogroup geometry geometry
///@{

/// \class colmesh
/// A colmesh, or "column mesh", is a semi-structured mesh consisting
/// of a set of columns of cells made by extruding polygons along the z axis.
/// Accordingly, each cell is a prism with polygonal "top" and "bottom" faces
/// whose normals align with the z axis, and rectangular "lateral" faces with
/// normals that lie in the xy plane.
/// To help with load balancing, a colmesh is divided into chunks. The chunks
/// consist of load-balanced sets of columns that span a well-defined segment
/// along the z axis.
typedef struct colmesh_t colmesh_t;

/// \enum colmesh_centering_t
/// Centerings for data on colmeshes. The structure of a colmesh
/// distinguishes between elements aligned with the z axis and those that
/// lie in the x-y plane.
typedef enum
{
  COLMESH_CELL = 0,
  COLMESH_XYFACE = 1,
  COLMESH_ZFACE = 2,
  COLMESH_XYEDGE = 3,
  COLMESH_ZEDGE = 4,
  COLMESH_NODE = 5
} colmesh_centering_t;

typedef struct colmesh_column_t colmesh_column_t;

/// \class colmesh_chunk
/// A group of columns within a colmesh.
struct colmesh_chunk_t
{
  /// The number of interior (polygonal) columns in this chunk.
  int num_columns;

  /// The number of ghost (polygonal) columns in this chunk.
  int num_ghost_columns;

  /// The number of vertical (z) cells in this chunk. Note that
  /// * the number of "z" faces is \ref num_zcells + 1.
  /// * the number of "z" edges is \ref num_zcells.
  /// * the number of "z" nodes is \ref num_zcells + 1.
  int num_z_cells;

  /// The z coordinate of the lower boundary of the chunk.
  real_t z1;

  /// The z coordinate of the upper boundary of the chunk.
  real_t z2;

  /// Offsets of lateral (xy) faces attached to columns, stored in compressed-row
  /// storage (CRS) format. column_xy_face_offsets[i] stores the offset within
  /// \ref column_faces for the ith xy face.
  int* column_xy_face_offsets;

  /// The indices of xy faces for columns, stored in CRS format.
  int* column_xy_faces;

  /// The columns attached to the xy faces in the mesh. Each xy face
  /// connects 2 columns, so the first column of the ith face is
  /// face_columns[2*i] and the second is face_columns[2*i+1].
  int* xy_face_columns;

  /// The total number of lateral (xy) faces at a single z location.
  int num_xy_faces;

  /// The total number of lateral (xy) edges at a single z location.
  int num_xy_edges;

  /// The nodes attached to an xy edge. The first node for an edge e is
  /// xy_edge_nodes[2*e], and the second is xy_edge_nodes[2*e+1].
  int* xy_edge_nodes;

  /// The total number of nodes at a single z location.
  int num_xy_nodes;

  /// The positions of the nodes in this chunk (in the xy plane).
  point2_t* xy_nodes;
};
typedef struct colmesh_chunk_t colmesh_chunk_t;

/// \class colmesh_fragment
/// This opaque type holds a planar_polymesh and send/receive maps that connect it
/// to other planar polymesh objects.
typedef struct colmesh_fragment_t colmesh_fragment_t;

/// \class colmesh_fragment_map
/// An distributed unordered map of xy indices to planar polymesh objects
/// representing fragments of a global planar polymesh. Used for constructing
/// especially large colmeshes.
DEFINE_UNORDERED_MAP(colmesh_fragment_map, int, colmesh_fragment_t*, int_hash, int_equals)

/// Adds a fragment of a planar polymesh to the given colmesh_fragment_map, along with send and
/// receive maps that define how the fragment connects to others.
/// \param [in] xy_index The xy index of the fragment within the global planar polymesh.
/// \param [in] fragment A planar polymesh representing the fragment itself. This argument is consumed
///                      by the map.
/// \param [in] send_map An exchanger_proc_map that connects send cells in this fragment to
///                      other fragments on this or other processes. Consumed by the map.
/// \param [in] receive_map An exchanger_proc_map that connects receive cells in this fragment to
///                         other fragments on this or other processes. Consumed by the map.
/// \memberof colmesh_fragment_map
void colmesh_fragment_map_add(colmesh_fragment_map_t* map,
                              int xy_index,
                              planar_polymesh_t* fragment,
                              exchanger_proc_map_t* send_map,
                              exchanger_proc_map_t* receive_map);

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct colmeshs.
// \ref colmesh_finalize must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Creates a new empty colmesh defined within the segment [z1, z2] extruded
/// from the given planar mesh of polygonal "column" cells.
/// \param [in] comm The communicator on which the mesh is constructed.
/// \param [in] columns A planar polygonal mesh that defines a set of connected
///                     polygonal columns for the colmesh. Must be the same mesh on all
///                     processes.
/// \param [in] z1 The z coordinate of the lower boundary of the mesh.
/// \param [in] z2 The z coordinate of the upper boundary of the mesh.
/// \param [in] num_xy_chunks The number of chunks in the distributed mesh within the xy plane.
/// \param [in] num_z_chunks The number of chunks in the distributed mesh along the z axis.
/// \param [in] nz_per_chunk The number of mesh cells along the z axis in a chunk.
/// \param [in] periodic_in_z True if the mesh is periodic along the z axis, false if not.
/// \returns A newly created colmesh containing no chunks.
/// \memberof colmesh
colmesh_t* create_empty_colmesh(MPI_Comm comm,
                                planar_polymesh_t* columns,
                                real_t z1, real_t z2,
                                int num_xy_chunks, int num_z_chunks,
                                int nz_per_chunk, bool periodic_in_z);

/// Creates a new empty colmesh defined within the segment [z1, z2] extruded
/// from the a planar mesh that forms the union of a distributed set of
/// fragments.
/// \param [in] comm The communicator on which the mesh is constructed.
/// \param [in] local_fragments An unordered map of integer xy indices to
///             planar_polymesh objects that represent locally-stored fragments
///             of the (global) planar mesh. This is a distributed map: each process
///             stores one or more fragments locally. The fragments on different
///             processes may overlap, and each one stores send/receive data that allows
///             the distributed global mesh to be constructed. This map is consumed
///             by this function.
/// \param [in] z1 The z coordinate of the lower boundary of the mesh.
/// \param [in] z2 The z coordinate of the upper boundary of the mesh.
/// \param [in] num_xy_chunks The total number of chunks within the xy plane.
/// \param [in] num_z_chunks The total number of chunks along the z axis.
/// \param [in] nz_per_chunk The total number of cells along the z axis for each chunk.
/// \param [in] periodic_in_z True if the mesh is periodic along the z axis, false if not.
/// \returns A newly created colmesh containing no chunks.
/// \memberof colmesh
/// \collective Collective on comm.
colmesh_t* create_empty_colmesh_from_fragments(MPI_Comm comm,
                                               colmesh_fragment_map_t* local_fragments,
                                               real_t z1, real_t z2,
                                               int num_xy_chunks,
                                               int num_z_chunks,
                                               int nz_per_chunk,
                                               bool periodic_in_z);

/// Inserts a new locally-stored chunk with the given xy and z indices into the mesh.
/// \param [in] xy_index The index identifying the polygonal column that contains the new chunk.
/// \param [in] z_index The index identifying the vertical (z) segment that contains the new chunk.
/// \memberof colmesh
void colmesh_insert_chunk(colmesh_t* mesh, int xy_index, int z_index);

/// Finalizes the construction process for the colmesh. This must be called
/// before any of the mesh's usage methods (below) are invoked. Should only
/// be called once.
/// \memberof colmesh
/// \collective Collective on the mesh's communicator.
void colmesh_finalize(colmesh_t* mesh);

//------------------------------------------------------------------------
//                          Ready-made constructors
//------------------------------------------------------------------------
// The following methods create colmeshes that are ready to go.
// No need to call colmesh_finalize() on these.
//------------------------------------------------------------------------

/// Creates a colmesh consisting of polygonal columns from the given
/// (global) planar polygonal mesh. The partitioning does not minimize
/// communication, so you might want to call \ref repartition_colmesh on the
/// resulting mesh.
/// \param [in] comm The communicator on which the mesh is constructed.
/// \param [in] columns A planar polygonal mesh that defines a set of connected
///                     polygonal columns for the colmesh. Consumed by this
///                     function.
/// \param [in] z1 The z coordinate of the lower boundary of the mesh.
/// \param [in] z2 The z coordinate of the upper boundary of the mesh.
/// \param [in] nz The total number of cells along the z axis.
/// \param [in] periodic_in_z True if the mesh is periodic along the z axis, false if not.
/// \returns A newly created colmesh.
/// \memberof colmesh
/// \collective Collective on comm.
colmesh_t* colmesh_new(MPI_Comm comm,
                       planar_polymesh_t* columns,
                       real_t z1, real_t z2,
                       int nz, bool periodic_in_z);

//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a colmesh has been
// fully constructed and finalized.
//------------------------------------------------------------------------

/// Returns true if the given mesh is finalized, false otherwise.
/// \memberof colmesh
bool colmesh_is_finalized(colmesh_t* mesh);

/// Destroys the given mesh.
/// \memberof colmesh
void colmesh_free(colmesh_t* mesh);

/// Retrieves metadata for chunks in this colmesh.
/// \param [out] num_xy_chunks The number of chunks within the xy plane.
/// \param [out] num_z_chunks The number of chunks along the z axis.
/// \param [out] nz_per_chunk The number of vertical cells per chunk.
/// \memberof colmesh
void colmesh_get_chunk_info(colmesh_t* mesh,
                            int* num_xy_chunks,
                            int* num_z_chunks,
                            int* nz_per_chunk);

/// Retrieves metadata for the z axis in this colmesh.
/// \param [out] z1 Stores the z coordinate of the lowest plane in the mesh.
/// \param [out] z2 Stores the z coordinate of the highest plane in the mesh.
/// \param [out] periodic Stores whether the z axis is periodic.
/// \memberof colmesh
void colmesh_get_z_info(colmesh_t* mesh,
                        real_t* z1,
                        real_t* z2,
                        bool* periodic);

/// Returns true if this colmesh is topologically correct, false if not.
/// \param [in] reason If non-NULL, stores a pointer to an internal string
///                    explaining why the mesh isn't topologically correct.
///                    Only used if the function returns false.
/// \memberof colmesh
bool colmesh_is_valid(colmesh_t* mesh, char** reason);

/// Returns the MPI communicator on which the colmesh is defined.
/// \memberof colmesh
MPI_Comm colmesh_comm(colmesh_t* mesh);

/// Returns the number of locally-stored chunks in the colmesh.
/// \memberof colmesh
int colmesh_num_chunks(colmesh_t* mesh);

/// Returns true if the mesh has a locally-stored chunk with the given
/// xy and z indices.
/// \param [in] xy_index The xy index of the chunk in question.
/// \param [in] z_index The z index of the chunk in question.
bool colmesh_has_chunk(colmesh_t* mesh, int xy_index, int z_index);

/// Returns the locally-stored chunk at the given xy and z indices, or NULL
/// this chunk is not locally-stored.
/// \param [in] xy_index The xy index of the chunk in question.
/// \param [in] z_index The z index of the chunk in question.
colmesh_chunk_t* colmesh_chunk(colmesh_t* mesh, int xy_index, int z_index);

/// Traverses the locally-stored chunks in the mesh.
/// \param [in,out] pos Controls the traversal. Set to 0 to reset traversal.
/// \param [out] xy_index Stores the xy index of the next chunk.
/// \param [out] z_index Stores the z index of the next chunk.
/// \param [out] chunk Stores the next chunk.
/// \returns true if more locally-stored chunks remain, false otherwise.
/// \memberof colmesh
bool colmesh_next_chunk(colmesh_t* mesh, int* pos,
                        int* xy_index, int* z_index,
                        colmesh_chunk_t** chunk);

/// Returns true if the given colmesh chunk is topologically correct, false if
/// not.
/// \param [out] reason If non-NULL, stores the reason that the colmesh chunk
///                     isn't topologically correct. Only used if this function
///                     returns false.
/// \memberof colmesh_chunk
bool colmesh_chunk_is_valid(colmesh_chunk_t* chunk, char** reason);

/// Returns a newly created polygon that represents the geometry of the
/// given column in the chunk.
/// \memberof colmesh_chunk
polygon_t* colmesh_chunk_polygon(colmesh_chunk_t* chunk, int column);

/// Returns the number of xy faces for the given column in the chunk.
/// \memberof colmesh_chunk
static inline int colmesh_chunk_column_num_xy_faces(colmesh_chunk_t* chunk,
                                                    int column)
{
  return chunk->column_xy_face_offsets[column+1] - chunk->column_xy_face_offsets[column];
}

/// Returns the indices of the xy faces for the given column in the chunk.
/// \param [in] column The index for the column.
/// \param [out] xy_faces An array big enough to store the (xy) indices of the
///                       xy faces of the column.
/// \memberof colmesh_chunk
static inline void colmesh_chunk_column_get_xy_faces(colmesh_chunk_t* chunk,
                                                     int column,
                                                     int* xy_faces)
{
  int start = chunk->column_xy_face_offsets[column];
  int end = chunk->column_xy_face_offsets[column+1];
  for (int f = start; f < end; ++f)
  {
    int face = chunk->column_xy_faces[f];
    if (face < 0)
      face = ~face;
    xy_faces[f-start] = face;
  }
}

/// Returns the "oriented" indices of the xy faces for the given column in the chunk.
/// Returns a non-negative index for xy faces whose edges/nodes are to be traversed in
/// order, and the one's complement for xy faces for which the traversal order is reversed.
/// \param [in] column The index for the column.
/// \param [out] xy_faces An array big enough to store the (xy) indices of the
///                       xy faces of the column.
/// \memberof colmesh_chunk
static inline void colmesh_chunk_column_get_oriented_xy_faces(colmesh_chunk_t* chunk,
                                                              int column,
                                                              int* xy_faces)
{
  int start = chunk->column_xy_face_offsets[column];
  int end = chunk->column_xy_face_offsets[column+1];
  for (int f = start; f < end; ++f)
    xy_faces[f-start] = chunk->column_xy_faces[f];
}

/// Traverses the neighboring columns of the given column in a chunk.
/// order, and the one's complement for xy faces for which the traversal order is reversed.
/// \param [in] column The index for the column.
/// \param [in,out] pos Controls the traversal. Set to 0 to reset traversal.
/// \param [out] neighboring_column Stores the index of the next neighboring column.
/// \returns true if more neighboring columns remain, false if the traversal is complete.
/// \memberof colmesh_chunk
static inline bool colmesh_chunk_column_next_neighbor(colmesh_chunk_t* chunk,
                                                      int column,
                                                      int* pos,
                                                      int* neighbor_column)
{
  int start = chunk->column_xy_face_offsets[column];
  int end = chunk->column_xy_face_offsets[column];
  if (*pos < (end - start))
  {
    int face = chunk->column_xy_faces[start+(*pos)];
    if (face < 0)
      *neighbor_column = chunk->xy_face_columns[2*(~face)+1];
    else
      *neighbor_column = chunk->xy_face_columns[2*face];
    ++(*pos);
    return true;
  }
  else
    return false;
}

/// Returns the xy and z indices for the nodes of the given xy face at the
/// given z index. The nodes of an xy face are ordered so that they're traversed
/// counterclockwise to give a right-handed orientation to the face.
/// Here's how we do it:
/// * Nodes 0 and 1 are the nodes of the top edge for the xy face,
///   traversed clockwise around the polygon forming the top Z face.
/// * Node 2 is the node on the far side of the z edge connected to node 1
///   for this xy face.
/// * Node 3 is the node on the far side of the xy edge connected to node 2
///   for this xy face. This should be the node connected to node 0 by a z
///   edge, attached to this face.
/// \param xy_face_index The index of the xy face within this chunk.
/// \param z_index The index indicating the z position of the xy face.
/// \param node_xy_indices An array of length 4 that stores the xy indices
///                        of the nodes for this face.
/// \param node_z_indices An array of length 4 that stores the z indices
///                       of the nodes for this face.
/// \memberof colmesh_chunk
static inline void colmesh_chunk_xy_face_get_nodes(colmesh_chunk_t* chunk,
                                                   int xy_face_index,
                                                   int z_index,
                                                   int* node_xy_indices,
                                                   int* node_z_indices)
{
  node_xy_indices[0] = chunk->xy_edge_nodes[2*xy_face_index];
  node_xy_indices[1] = chunk->xy_edge_nodes[2*xy_face_index+1];
  node_xy_indices[2] = chunk->xy_edge_nodes[2*xy_face_index+1];
  node_xy_indices[3] = chunk->xy_edge_nodes[2*xy_face_index];
  node_z_indices[0] = z_index+1;
  node_z_indices[1] = z_index+1;
  node_z_indices[2] = z_index;
  node_z_indices[3] = z_index;
}

/// Returns the xy and z indices for the edges of the given xy face at the
/// given z index. The edges of an xy face are ordered so that they're traversed
/// counterclockwise to give a right-handed orientation to the face.
/// Here's how we do it:
/// * Edge 0 is the top edge for the xy face,
/// * Edge 1 is the edge connecting the face's top edge to its bottom edge,
///   traversed counterclockwise from edge 0.
/// * Edge 2 is the bottom edge for the xy face,
/// * Edge 3 is the remaining edge for the face, opposite edge 1.
/// \param xy_face_index The index of the xy face within this chunk.
/// \param z_index The index indicating the z position of the xy face.
/// \param edge_xy_indices An array of length 4 that stores the xy indices
///                        of the edges for this face.
/// \param edge_z_indices An array of length 4 that stores the z indices
///                       of the edges for this face.
/// \memberof colmesh_chunk
static inline void colmesh_chunk_xy_face_get_edges(colmesh_chunk_t* chunk,
                                                   int xy_face_index,
                                                   int z_index,
                                                   int* edge_xy_indices,
                                                   int* edge_z_indices)
{
  edge_xy_indices[0] = xy_face_index;
  edge_xy_indices[1] = xy_face_index+1;
  edge_xy_indices[2] = xy_face_index+1;
  edge_xy_indices[3] = xy_face_index;
  edge_z_indices[0] = z_index+1;
  edge_z_indices[1] = z_index+1;
  edge_z_indices[2] = z_index;
  edge_z_indices[3] = z_index;
}

/// Returns the number of nodes for the given z face in the chunk.
/// \memberof colmesh_chunk
static inline int colmesh_chunk_z_face_num_nodes(colmesh_chunk_t* chunk,
                                                 int z_face)
{
  // There are as many nodes on a z face as there are xy faces on a column.
  return colmesh_chunk_column_num_xy_faces(chunk, z_face);
}

/// Returns the indices of the nodes for the given z face in the chunk.
/// \param [in] z_face The index for the z face (same as the column index).
/// \param [out] nodes An array big enough to store the indices of the
///                    nodes of the z face.
/// \memberof colmesh_chunk
void colmesh_chunk_z_face_get_nodes(colmesh_chunk_t* chunk,
                                    int z_face,
                                    int* nodes);

/// Retrieves the coordinates of the node with the given xy and z indices in this
/// chunk.
/// \param [in] xy_index The index identifying the position of the node in the xy plane.
/// \param [in] z_index The index identifying the position of the node along the z axis.
/// \param [out] node_pos Stores the coordinates of the node.
void colmesh_chunk_get_node(colmesh_chunk_t* chunk,
                            int xy_index, int z_index,
                            point_t* node_pos);

typedef struct colmesh_field_t colmesh_field_t;

/// Repartitions the given colmesh and redistributes data to each of the
/// given fields. Here, the old meshes and fields are consumed, and new ones
/// are created in their place. Weights can be provided for each patch, and
/// the partitioning is performed so that the load imbalance does not exceed
/// the given tolerance.
/// \note In addition, each repartitioned field needs to have any boundary
/// conditions reinstated, since these boundary conditions are not
/// transmitted between processes.
/// \relates colmesh
void repartition_colmesh(colmesh_t** mesh,
                         int* weights,
                         real_t imbalance_tol,
                         colmesh_field_t** fields,
                         size_t num_fields);

typedef struct polymesh_t polymesh_t;

///@}

#endif

