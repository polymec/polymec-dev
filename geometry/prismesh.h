// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRISMESH_H
#define POLYMEC_PRISMESH_H

#include "core/point.h"
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
  prismesh_CELL = 0,
  prismesh_XYFACE = 1,
  prismesh_ZFACE = 2,
  prismesh_XYEDGE = 3,
  prismesh_ZEDGE = 3,
  prismesh_NODE = 3
} prismesh_centering_t;

typedef struct prismesh_column_t prismesh_column_t;

/// \class prismesh_column
/// A column of cells within a prismesh.
struct prismesh_column_t
{
  /// Index of the column within the prismesh.
  size_t index;
  /// Polygon representing the top and bottom facet.
  polygon_t* polygon;
  /// Number of cells in the column.
  size_t num_cells; 
  /// Array of neighboring columns (arranged in the same order as the 
  /// vertices/edges of the polygon). NULL if the column is on a parallel 
  /// domain boundary or on the boundary of the mesh.
  prismesh_column_t** neighbors;
}; 

DEFINE_ARRAY(prism_column_array, prismesh_column_t*)

/// \class prismesh_layer
/// A layer of columns within a prismesh.
struct prismesh_layer_t 
{
  /// The underlying mesh.
  prismesh_t* mesh;

  /// The z coordinate of the layer's lower boundary.
  real_t z1;
  /// The z coordinate of the layer's upper boundary.
  real_t z2;

  /// Columns locally stored within this layer.
  prism_column_array_t* columns;

  /// cell_xyfaces[m][n] stores the nth xy-face for the mth cell.
  int** cell_xyfaces;

  /// xyface_zedges[m][n] stores the nth z-edge for the mth xy-face.
  int** xyface_zedges;

  /// zface_xyedges[m][n] stores the nth xy-edge for the mth z-face.
  int** zface_xyedges;

  /// xyface_nodes[m][n] stores the nth node for the mth xy-face.
  int** xyface_nodes;

  /// zface_nodes[m][n] stores the nth node for the mth z-face.
  int** zface_nodes;
};
typedef struct prismesh_layer_t prismesh_layer_t;

/// Creates a prismesh with the given number of layers and columns. 
/// Layers and columns are automatically created and are managed by the mesh.
/// \param comm [in] The MPI communicator used by the mesh.
/// \param num_columns [in] The number of columns in the mesh.
/// \param num_vertical_cells [in] The number of cells along the z axis.
/// \param z [in] An array of length (num_vertical+1) that gives the z 
///               coordinates of the interfaces between cells, from bottom
///               to top.
/// \memberof prismesh
prismesh_t* prismesh_new(MPI_Comm comm, 
                         size_t num_columns, 
                         size_t num_vertical_cells,
                         real_t* z);
 
//------------------------------------------------------------------------
//                        Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct prismeshes.
// prismesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Sets the geometry and the neighbors for the given column in the mesh. 
/// \memberof prismesh
/// \param [in] column The index for the desired column.
/// \param [in] polygon The polygonal face for the top and bottom of each 
///                     cell in the column.
/// \param [in] neighbors An array of column indices indicating the neighbors 
///                       of the given column. The length of this array is 
///                       \ref polygon_num_edges(polygon).
void prismesh_set_column(prismesh_t* mesh, 
                         size_t column, 
                         polygon_t* polygon,
                         size_t* neighbors);

/// Finalizes the construction process for the mesh. This must be called 
/// before any of the mesh's usage methods (below) are invoked. Should only 
/// be called once.
/// \memberof prismesh
void prismesh_finalize(prismesh_t* mesh);

//------------------------------------------------------------------------
//                        Usage methods
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

/// Returns the number of layers in the prismesh.
/// \memberof prismesh
size_t prismesh_num_layers(prismesh_t* mesh);

/// Returns the total number of columns in the prismesh, as seen 
/// looking down on the top layer. This is *not* the same as the 
/// total number of prismesh_columns globally accessible in the mesh, 
/// since these columns are divided along the z axis for scalability.
/// \memberof prismesh
size_t prismesh_num_columns(prismesh_t* mesh);

/// Returns the total number of vertical cells in each column in the prismesh.
/// \memberof prismesh
size_t prismesh_num_vertical_cells(prismesh_t* mesh);

/// Returns the polygon associated with the given column.
/// \memberof prismesh
polygon_t* prismesh_polygon(prismesh_t* mesh, size_t column);

/// Traverses the locally-stored layers in the mesh, returning true and the 
/// next layer if the traversal is incomplete, false otherwise. 
/// \memberof prismesh
bool prismesh_next_layer(prismesh_t* mesh, int* pos, prismesh_layer_t** layer);

/// Returns the number of columns in the prismesh layer.
/// \memberof prismesh_layer
size_t prismesh_layer_num_columns(prismesh_layer_t* layer);

/// Returns the bottom and top z coordinates of the prismesh layer.
/// \memberof prismesh_layer
/// \param [out] z1 Stores the bottom z coordinate of the layer.
/// \param [out] z2 Stores the top z coordinate of the layer.
void prismesh_layer_get_bounds(prismesh_layer_t* layer, real_t* z1, real_t *z2);

/// Traverses the locally-stored columns in the prismesh layer, returning true 
/// and the next column if the traversal is incomplete, false otherwise. 
/// \memberof prismesh_layer
bool prismesh_layer_next_column(prismesh_layer_t* layer, int* pos, prismesh_column_t** column);

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

