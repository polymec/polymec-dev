// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRISMESH_H
#define POLYMEC_PRISMESH_H

#include "core/point.h"
#include "geometry/planar_polymesh.h"
#include "geometry/polygon.h"

/// \addtogroup geometry geometry
///@{

/// \class prismesh
/// A prismesh, or polygonal extruded mesh, is a semi-structured mesh consisting 
/// of a set of columns of cells. Each cell has "top" and "bottom" faces that 
/// are polygons, and rectangular lateral faces. The "top" and "bottom" faces 
/// lie in the x-y plane, and the column extends along the z axis. The columns 
/// are divided into layers along the z axis to facilitate load balancing of 
/// large prismeshes.
typedef struct prismesh_t prismesh_t;

/// \enum prismesh_centering
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

/// \class prismesh_layer
/// A layer of columns within a prismesh.
struct prismesh_layer_t 
{
  /// The number of (polygonal) columns in this layer.
  size_t num_columns;

  /// The number of vertical (z) cells in this layer. Note that 
  /// * the number of "z" faces is \ref num_zcells + 1.
  /// * the number of "z" edges is \ref num_zcells.
  /// * the number of "z" nodes is \ref num_zcells + 1.
  size_t num_z_cells;

  /// The z coordinate of the lower boundary of the layer.
  real_t z1;

  /// The z coordinate of the upper boundary of the layer.
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

  /// The positions of the nodes in this layer (in the xy plane).
  point2_t* xy_nodes;
};
typedef struct prismesh_layer_t prismesh_layer_t;

/// Creates a prismesh consisting of polygonal columns from the given 
/// planar polygonal mesh on the same communicator as that mesh. The columns 
/// are distributed in the same way as the planar polygonal mesh, which is 
/// probably not optimal for calculations on the prismesh.
/// \param columns [in] A planar polygonal mesh that defines a set of connected
///                     polygonal columns for the prismesh.
/// \param num_vertical_cells [in] The number of cells along the z axis.
/// \param z1 [in] The z coordinate of the lower boundary of the mesh.
/// \param z2 [in] The z coordinate of the upper boundary of the mesh.
/// \returns A newly created prismesh.
/// \memberof prismesh
prismesh_t* prismesh_new(planar_polymesh_t* columns,
                         size_t num_vertical_cells,
                         real_t z1, real_t z2);
 
/// Destroys the given mesh.
/// \memberof prismesh
void prismesh_free(prismesh_t* mesh);

/// Returns the MPI communicator on which the prismesh is defined.
/// \memberof prismesh
MPI_Comm prismesh_comm(prismesh_t* mesh);

/// Returns the number of layers in the prismesh.
/// \memberof prismesh
size_t prismesh_num_layers(prismesh_t* mesh);

/// Returns the total number of columns in the prismesh, as seen 
/// looking down on the top layer. This is *not* the same as the 
/// total number of prismesh_columns globally accessible in the mesh, 
/// since these columns are divided along the z axis for scalability.
/// \memberof prismesh
size_t prismesh_num_columns(prismesh_t* mesh);

/// Returns the total number of cells along the z axis in the prismesh.
/// \memberof prismesh
size_t prismesh_num_vertical_cells(prismesh_t* mesh);

/// Returns the total number of locally-stored cells in the prismesh.
/// \memberof prismesh
size_t prismesh_num_cells(prismesh_t* mesh);

/// Returns the z coordinate of the bottom boundary of the mesh.
real_t primesh_z1(prismesh_t* mesh);

/// Returns the z coordinate of the top boundary of the mesh.
real_t primesh_z2(prismesh_t* mesh);

/// Returns a newly created polygon that represents the geometry of the 
/// given column.
/// \memberof prismesh
polygon_t* prismesh_polygon(prismesh_t* mesh, size_t column);

/// Traverses the locally-stored layers in the mesh, returning true and the 
/// next layer if the traversal is incomplete, false otherwise. 
/// \memberof prismesh
bool prismesh_next_layer(prismesh_t* mesh, int* pos, prismesh_layer_t** layer);

/// Returns the number of xy faces for the given column in the layer.
/// \memberof prismesh_layer
static inline int prismesh_layer_column_num_xy_faces(prismesh_layer_t* layer,
                                                     int column)
{
  return layer->column_xy_face_offsets[column+1] - layer->column_xy_face_offsets[column];
}

/// Returns the indices of the xy faces for the given column in the layer.
/// \param column [in] The index for the column.
/// \param xy_faces [out] An array big enough to store the (xy) indices of the
///                       xy faces of the column.
/// \memberof prismesh_layer
static inline void prismesh_layer_column_get_xy_faces(prismesh_layer_t* layer,
                                                      int column,
                                                      int* xy_faces)
{
  int start = layer->column_xy_face_offsets[column];
  int end = layer->column_xy_face_offsets[column+1];
  for (int f = start; f < end; ++f)
    xy_faces[f] = layer->column_xy_faces[f];
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
/// \param xy_face_index The index of the xy face within this layer.
/// \param z_index The index indicating the z position of the xy face.
/// \param node_xy_indices An array of length 4 that stores the xy indices
///                        of the nodes for this face.
/// \param node_z_indices An array of length 4 that stores the z indices
///                       of the nodes for this face.
/// \memberof prismesh_layer
static inline void prismesh_layer_xy_face_get_nodes(prismesh_layer_t* layer,
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
/// \param xy_face_index The index of the xy face within this layer.
/// \param z_index The index indicating the z position of the xy face.
/// \param edge_xy_indices An array of length 4 that stores the xy indices
///                        of the edges for this face.
/// \param edge_z_indices An array of length 4 that stores the z indices
///                       of the edges for this face.
/// \memberof prismesh_layer
static inline void prismesh_layer_xy_face_get_edges(prismesh_layer_t* layer,
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

