// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PEXMESH_H
#define POLYMEC_PEXMESH_H

#include "core/point.h"
#include "geometry/polygon.h"

/// \addtogroup geometry geometry
///@{

/// \class pexmesh
/// A pexmesh, or polygonal extruded mesh, is a semi-structured mesh consisting 
/// of a set of columns of cells. Each cell has "top" and "bottom" faces that 
/// are polygons, and rectangular lateral faces. The "top" and "bottom" faces 
/// lie in the x-y plane, and the column extends along the z axis. The columns 
/// are divided into layers along the z axis to facilitate load balancing of 
/// large pexmeshes.
typedef struct pexmesh_t pexmesh_t;

/// \enum pexmesh_centering
/// Centerings for data on pexmeshes.
typedef enum
{
  PEXMESH_CELL = 0,
  PEXMESH_FACE = 1,
  PEXMESH_EDGE = 2,
  PEXMESH_NODE = 3
} pexmesh_centering_t;

/// \class pexmesh_layer
/// A layer of columns within a pexmesh.
typedef struct pexmesh_layer_t pexmesh_layer_t;

typedef struct pexmesh_column_t pexmesh_column_t;

/// \class pexmesh_column
/// A column of cells within a pexmesh.
struct pexmesh_column_t
{
  /// Index of the column within the pexmesh.
  size_t index;
  /// Polygon representing the top and bottom facet.
  polygon_t* polygon;
  /// Number of cells in the column.
  size_t num_cells; 
  /// Array of neighboring columns (arranged in the same order as the 
  /// vertices/edges of the polygon). NULL if the column is on a parallel 
  /// domain boundary or on the boundary of the mesh.
  pexmesh_column_t** neighbors;
}; 

/// Creates a pexmesh with the given number of layers and columns. 
/// Layers and columns are automatically created and are managed by the mesh.
/// \param comm [in] The MPI communicator used by the mesh.
/// \param num_columns [in] The number of columns in the mesh.
/// \param num_vertical_cells [in] The number of cells along the z axis.
/// \param z [in] An array of length (num_vertical+1) that gives the z 
///               coordinates of the interfaces between cells, from bottom
///               to top.
/// \memberof pexmesh
pexmesh_t* pexmesh_new(MPI_Comm comm, 
                       size_t num_columns, 
                       size_t num_vertical_cells,
                       real_t* z);

//------------------------------------------------------------------------
//                        Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct pexmeshes.
// pexmesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Sets the geometry and the neighbors for the given column in the mesh. 
/// \memberof pexmesh
/// \param [in] column The index for the desired column.
/// \param [in] polygon The polygonal face for the top and bottom of each 
///                     cell in the column.
/// \param [in] neighbors An array of column indices indicating the neighbors 
///                       of the given column. The length of this array is 
///                       \ref polygon_num_edges(polygon).
void pexmesh_set_column(pexmesh_t* mesh, 
                        size_t column, 
                        polygon_t* polygon,
                        size_t* neighbors);

/// Finalizes the construction process for the mesh. This must be called 
/// before any of the mesh's usage methods (below) are invoked. Should only 
/// be called once.
/// \memberof pexmesh
void pexmesh_finalize(pexmesh_t* mesh);

//------------------------------------------------------------------------
//                        Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a pexmesh has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

/// Destroys the given mesh.
/// \memberof pexmesh
void pexmesh_free(pexmesh_t* mesh);

/// Returns the MPI communicator on which the pexmesh is defined.
/// \memberof pexmesh
MPI_Comm pexmesh_comm(pexmesh_t* mesh);

/// Returns the number of layers in the pexmesh.
/// \memberof pexmesh
size_t pexmesh_num_layers(pexmesh_t* mesh);

/// Returns the total number of columns in the pexmesh, as seen 
/// looking down on the top layer. This is *not* the same as the 
/// total number of pexmesh_columns globally accessible in the mesh, 
/// since these columns are divided along the z axis for scalability.
/// \memberof pexmesh
size_t pexmesh_num_columns(pexmesh_t* mesh);

/// Returns the total number of vertical cells in each column in the pexmesh.
/// \memberof pexmesh
size_t pexmesh_num_vertical_cells(pexmesh_t* mesh);

/// Returns the polygon associated with the given column.
/// \memberof pexmesh
polygon_t* pexmesh_polygon(pexmesh_t* mesh, size_t column);

/// Traverses the locally-stored layers in the mesh, returning true and the 
/// next layer if the traversal is incomplete, false otherwise. 
/// \memberof pexmesh
bool pexmesh_next_layer(pexmesh_t* mesh, int* pos, pexmesh_layer_t** layer);

/// Returns the number of columns in the pexmesh layer.
/// \memberof pexmesh_layer
size_t pexmesh_layer_num_columns(pexmesh_layer_t* layer);

/// Returns the bottom and top z coordinates of the pexmesh layer.
/// \memberof pexmesh_layer
/// \param [out] z1 Stores the bottom z coordinate of the layer.
/// \param [out] z2 Stores the top z coordinate of the layer.
void pexmesh_layer_get_bounds(pexmesh_layer_t* layer, real_t* z1, real_t *z2);

/// Traverses the locally-stored columns in the pexmesh layer, returning true 
/// and the next column if the traversal is incomplete, false otherwise. 
/// \memberof pexmesh_layer
bool pexmesh_layer_next_column(pexmesh_layer_t* layer, int* pos, pexmesh_column_t** column);

typedef struct pexmesh_field_t pexmesh_field_t;

/// Repartitions the given pexmesh and redistributes data to each of the 
/// given fields. Here, the old meshes and fields are consumed, and new ones 
/// are created in their place. Weights can be provided for each patch, and 
/// the partitioning is performed so that the load imbalance does not exceed 
/// the given tolerance.
/// \note In addition, each repartitioned field needs to have any boundary 
/// conditions reinstated, since these boundary conditions are not 
/// transmitted between processes.
/// \relates pexmesh
void repartition_pexmesh(pexmesh_t** mesh, 
                         int* weights,
                         real_t imbalance_tol,
                         pexmesh_field_t** fields,
                         size_t num_fields);

typedef struct polymesh_t polymesh_t;

/// Returns a polymesh that represents the same geometry as this pexmesh.
/// \memberof pexmesh
polymesh_t* pexmesh_as_polymesh(pexmesh_t* mesh);

///@}

#endif
