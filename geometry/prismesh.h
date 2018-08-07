// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRISMESH_H
#define POLYMEC_PRISMESH_H

#include "core/point.h"
#include "core/adj_graph.h"
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
  /// The index of the layer within this mesh.
  int index;

  /// The number of (polygonal) columns in this layer.
  size_t num_columns;

  /// The number of vertical (z) cells in this layer. Notice that 
  /// * the number of "z" faces is \ref num_zcells + 1.
  /// * the number of "z" edges is \ref num_zcells.
  /// * the number of "z" nodes is \ref num_zcells + 1.
  size_t num_z_cells;

  /// The z coordinate of the lower boundary of the layer.
  real_t z1;

  /// The z coordinate of the upper boundary of the layer.
  real_t z2;

  /// Offsets of lateral (xy-) faces attached to columns, stored in compressed-row 
  /// storage (CRS) format. column_xy_face_offsets[i] stores the offset within
  /// \ref column_faces for the ith xy-face.
  int* column_xy_face_offsets;

  /// The indices of xy-faces for columns, stored in CRS format.
  int* column_xy_faces;

  /// The columns attached to the xy-faces in the mesh. Each xy-face 
  /// connects 2 columns, so the first column of the ith face is 
  /// face_columns[2*i] and the second is face_columns[2*i+1].
  int* xy_face_columns;

  /// The total number of lateral (xy-) faces at a single z location.
  size_t num_xy_faces;

  /// The offsets of the sets of column nodes (polygon vertices) attached to 
  /// xy-faces, stored in CRS format.
  int* xy_face_node_offsets;

  /// The indices of column nodes (polygon vertices) attached to xy-faces, 
  /// stored in CRS format.
  int* xy_face_nodes;

  /// The offsets of the sets of column (polygon) edges attached to xy-faces, 
  /// stored in CRS format.
  int* xy_face_edge_offsets;
  /// The indices of column (polygon) xy-edges attached to xy-faces, stored in 
  /// CRS format.
  int* xy_face_edges;

  /// The total number of lateral (xy-) edges at a single z location.
  size_t num_xy_edges;

  /// The total number of nodes at a single z location.
  size_t num_xy_nodes;

  /// The positions of the nodes in this layer.
  point_t* nodes;
};
typedef struct prismesh_layer_t prismesh_layer_t;

/// Creates a prismesh consisting of the given polygonal columns, with
/// connectivity between columns defined by the given adjacency graph.
/// If the adjacency graph is defined on the MPI_COMM_SELF communicator, 
/// the prismesh is defined entirely on each process. Otherwise (if the 
/// adjacency graph is distributed across the processes for another 
/// communicator, the resulting prismesh is distributed in the same way.
/// In any case, the mesh emerging from this function is not necessarily 
/// optimally load balanced.
/// \param columns [in] An array of polygonal prism columns that define the 
///                     mesh.
/// \param num_columns [in] The length of the \ref columns array.
/// \param connectivity [in] An adjacency graph defining the connectivity 
///                          between columns. The vertices of this graph 
///                          are the polygonal columns, and the edges represent
///                          the faces that connect different columns.
/// \param num_vertical_cells [in] The number of cells along the z axis.
/// \param z1 [in] The z coordinate of the lower boundary of the mesh.
/// \param z2 [in] The z coordinate of the upper boundary of the mesh.
/// \memberof prismesh
prismesh_t* prismesh_new(polygon_t** columns,
                         size_t num_columns, 
                         adj_graph_t* connectivity,
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

/// Returns the polygon associated with the given column.
/// \memberof prismesh
polygon_t* prismesh_polygon(prismesh_t* mesh, size_t column);

/// Traverses the locally-stored layers in the mesh, returning true and the 
/// next layer if the traversal is incomplete, false otherwise. 
/// \memberof prismesh
bool prismesh_next_layer(prismesh_t* mesh, int* pos, prismesh_layer_t** layer);

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

