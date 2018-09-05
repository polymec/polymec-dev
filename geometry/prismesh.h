// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRISMESH_H
#define POLYMEC_PRISMESH_H

#include "core/point.h"
#include "core/exchanger.h"
#include "geometry/planar_polymesh.h"
#include "geometry/polygon.h"

/// \addtogroup geometry geometry
///@{

/// \class prismesh
/// A prismesh, or "prism mesh", is a semi-structured mesh consisting 
/// of a set of columns of cells made by extruding polygons along the z axis.
/// Accordingly, each cell has polygonal "top" and "bottom" faces whose normals 
/// align with the z axis, and rectangular "lateral" faces with normals that 
/// lie in the xy plane. 
/// To help with load balancing, a prismesh is divided into chunks. The chunks
/// consist of load-balanced sets of columns that span a well-defined segment 
/// along the z axis.
typedef struct prismesh_t prismesh_t;

/// \enum prismesh_centering_t
/// Centerings for data on prismeshes. The structure of a prismesh 
/// distinguishes between elements aligned with the z axis and those that 
/// lie in the x-y plane.
typedef enum
{
  PRISMESH_CELL = 0,
  PRISMESH_XYFACE = 1,
  PRISMESH_ZFACE = 2,
  PRISMESH_XYEDGE = 3,
  PRISMESH_ZEDGE = 3,
  PRISMESH_NODE = 3
} prismesh_centering_t;

typedef struct prismesh_column_t prismesh_column_t;

/// \class prismesh_chunk
/// A group of columns within a prismesh.
struct prismesh_chunk_t 
{
  /// The number of (polygonal) columns in this chunk.
  size_t num_columns;

  /// The number of vertical (z) cells in this chunk. Note that 
  /// * the number of "z" faces is \ref num_zcells + 1.
  /// * the number of "z" edges is \ref num_zcells.
  /// * the number of "z" nodes is \ref num_zcells + 1.
  size_t num_z_cells;

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
  size_t num_xy_faces;

  /// The total number of lateral (xy) edges at a single z location.
  size_t num_xy_edges;

  /// The total number of nodes at a single z location.
  size_t num_xy_nodes;

  /// The positions of the nodes in this chunk (in the xy plane).
  point2_t* xy_nodes;
};
typedef struct prismesh_chunk_t prismesh_chunk_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct prismeshs.
// prismesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Creates a new empty prismesh level defined within the segment [z1, z2] of 
/// of the z axis, with cellular columns defined by the given planar polymesh.
/// \param comm [in] The communicator on which the mesh is constructed.
/// \param columns [in] A planar polygonal mesh that defines a set of connected
///                     polygonal columns for the prismesh. Consumed by this 
///                     function.
/// \param z1 [in] The z coordinate of the lower boundary of the mesh.
/// \param z2 [in] The z coordinate of the upper boundary of the mesh.
/// \param num_xy_chunks [in] The number of chunks in the distributed mesh within the xy plane.
/// \param num_z_chunks [in] The number of chunks in the distributed mesh along the z axis.
/// \param nz_per_chunk [in] The number of mesh cells along the z axis in a chunk.
/// \returns A newly created prismesh containing no chunks.
/// \memberof prismesh
prismesh_t* create_empty_prismesh(MPI_Comm comm, 
                                  planar_polymesh_t* columns,
                                  real_t z1, real_t z2,
                                  size_t num_xy_chunks, size_t num_z_chunks,
                                  size_t nz_per_chunk);

/// Inserts a new locally-stored chunk with the given xy and z indices into the mesh.
/// \param xy_index [in] The index identifying the polygonal column that contains the new chunk.
/// \param z_index [in] The index identifying the vertical (z) segment that contains the new chunk.
/// \memberof prismesh
void prismesh_insert_chunk(prismesh_t* mesh, int xy_index, int z_index);

/// Finalizes the construction process for the prismesh. This must be called 
/// before any of the mesh's usage methods (below) are invoked. Should only 
/// be called once.
/// \memberof prismesh
void prismesh_finalize(prismesh_t* mesh);

//------------------------------------------------------------------------
//                          Ready-made constructors
//------------------------------------------------------------------------
// The following methods create prismeshes that are ready to go.
// No need to call prismesh_finalize() on these.
//------------------------------------------------------------------------

/// Creates a prismesh consisting of polygonal columns from the given 
/// (global) planar polygonal mesh. The partitioning does not minimize 
/// communication, so you might want to call repartition_prismesh on the 
/// resulting mesh.
/// \param comm [in] The communicator on which the mesh is constructed.
/// \param columns [in] A planar polygonal mesh that defines a set of connected
///                     polygonal columns for the prismesh. Consumed by this
///                     function.
/// \param z1 [in] The z coordinate of the lower boundary of the mesh.
/// \param z2 [in] The z coordinate of the upper boundary of the mesh.
/// \param nz [in] The total number of cells along the z axis.
/// \returns A newly created prismesh.
/// \memberof prismesh
/// \collective Collective on comm.
prismesh_t* prismesh_new(MPI_Comm comm,
                         planar_polymesh_t* columns,
                         real_t z1, real_t z2,
                         size_t nz);
 
//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a prismesh has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

/// Destroys the given mesh.
/// \memberof prismesh
void prismesh_free(prismesh_t* mesh);

/// Returns the MPI communicator on which the prismesh is defined.
/// \memberof prismesh
MPI_Comm prismesh_comm(prismesh_t* mesh);

/// Returns the number of locally-stored chunks in the prismesh.
/// \memberof prismesh
size_t prismesh_num_chunks(prismesh_t* mesh);

/// Returns the total number of chunks that the mesh can store in the xy plane.
/// \memberof prismesh
size_t prismesh_num_xy_chunks(prismesh_t* mesh);

/// Returns the total number of chunks that the mesh can store along the z axis.
/// \memberof prismesh
size_t prismesh_num_z_chunks(prismesh_t* mesh);

/// Returns the z coordinate of the bottom boundary of the mesh.
/// \memberof prismesh
real_t prismesh_z1(prismesh_t* mesh);

/// Returns the z coordinate of the top boundary of the mesh.
/// \memberof prismesh
real_t prismesh_z2(prismesh_t* mesh);

/// Returns true if the mesh has a locally-stored chunk with the given 
/// xy and z indices.
/// \param xy_index [in] The xy index of the chunk in question.
/// \param z_index [in] The z index of the chunk in question.
bool prismesh_has_chunk(prismesh_t* mesh, int xy_index, int z_index);

/// Returns the locally-stored chunk at the given xy and z indices, or NULL
/// this chunk is not locally-stored.
/// \param xy_index [in] The xy index of the chunk in question.
/// \param z_index [in] The z index of the chunk in question.
prismesh_chunk_t* prismesh_chunk(prismesh_t* mesh, int xy_index, int z_index);

/// Traverses the locally-stored chunks in the mesh.
/// \param pos [in,out] Controls the traversal. Set to 0 to reset traversal.
/// \param xy_index [out] Stores the xy index of the next chunk.
/// \param z_index [out] Stores the z index of the next chunk.
/// \param chunk [out] Stores the next chunk.
/// \returns true if more locally-stored chunks remain, false otherwise. 
/// \memberof prismesh
bool prismesh_next_chunk(prismesh_t* mesh, int* pos, 
                         int* xy_index, int* z_index,
                         prismesh_chunk_t** chunk);

/// Returns a newly created polygon that represents the geometry of the 
/// given column in the chunk.
/// \memberof prismesh_chunk
polygon_t* prismesh_chunk_polygon(prismesh_chunk_t* chunk, int column);

/// Returns the number of xy faces for the given column in the chunk.
/// \memberof prismesh_chunk
static inline int prismesh_chunk_column_num_xy_faces(prismesh_chunk_t* chunk,
                                                     int column)
{
  return chunk->column_xy_face_offsets[column+1] - chunk->column_xy_face_offsets[column];
}

/// Returns the indices of the xy faces for the given column in the chunk.
/// \param column [in] The index for the column.
/// \param xy_faces [out] An array big enough to store the (xy) indices of the
///                       xy faces of the column.
/// \memberof prismesh_chunk
static inline void prismesh_chunk_column_get_xy_faces(prismesh_chunk_t* chunk,
                                                      int column,
                                                      int* xy_faces)
{
  int start = chunk->column_xy_face_offsets[column];
  int end = chunk->column_xy_face_offsets[column+1];
  for (int f = start; f < end; ++f)
    xy_faces[f] = chunk->column_xy_faces[f];
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
/// \memberof prismesh_chunk
static inline void prismesh_chunk_xy_face_get_nodes(prismesh_chunk_t* chunk,
                                                    int xy_face_index,
                                                    int z_index,
                                                    int* node_xy_indices,
                                                    int* node_z_indices)
{
  node_xy_indices[0] = xy_face_index;
  node_xy_indices[1] = xy_face_index+1;
  node_xy_indices[2] = xy_face_index+1;
  node_xy_indices[3] = xy_face_index;
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
/// \memberof prismesh_chunk
static inline void prismesh_chunk_xy_face_get_edges(prismesh_chunk_t* chunk,
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
/// \memberof prismesh_chunk
static inline int prismesh_chunk_z_face_num_nodes(prismesh_chunk_t* chunk,
                                                  int z_face)
{
  // Nodes on z faces are indexed the same way as xy faces on columns.
  return prismesh_chunk_column_num_xy_faces(chunk, z_face);
}

/// Returns the indices of the nodes for the given z face in the chunk.
/// \param z_face [in] The index for the z face (same as the column index).
/// \param nodes [out] An array big enough to store the indices of the
///                    nodes of the z face.
/// \memberof prismesh_chunk
static inline void prismesh_chunk_z_face_get_nodes(prismesh_chunk_t* chunk,
                                                   int z_face,
                                                   int* nodes)
{
  // Nodes on z faces are indexed the same way as xy faces on columns.
  prismesh_chunk_column_get_xy_faces(chunk, z_face, nodes);
}

typedef struct prismesh_field_t prismesh_field_t;

/// Repartitions the given prismesh and redistributes data to each of the 
/// given fields. Here, the old meshes and fields are consumed, and new ones 
/// are created in their place. Weights can be provided for each patch, and 
/// the partitioning is performed so that the load imbalance does not exceed 
/// the given tolerance.
/// \note In addition, each repartitioned field needs to have any boundary 
/// conditions reinstated, since these boundary conditions are not 
/// transmitted between processes.
/// \relates prismesh
void repartition_prismesh(prismesh_t** mesh, 
                          int* weights,
                          real_t imbalance_tol,
                          prismesh_field_t** fields,
                          size_t num_fields);

typedef struct polymesh_t polymesh_t;

/// Returns a polymesh that represents the same geometry as this prismesh.
/// \memberof prismesh
polymesh_t* prismesh_as_polymesh(prismesh_t* mesh);

///@}

#endif

