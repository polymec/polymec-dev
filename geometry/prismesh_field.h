// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRISMESH_FIELD_H
#define POLYMEC_PRISMESH_FIELD_H

#include "geometry/prismesh.h"

/// \addtogroup geometry geometry
///@{

/// \class prismesh_field
/// This type represents field values defined on a prismesh.
typedef struct prismesh_field_t prismesh_field_t;

/// \struct prismesh_layer_data
/// This type holds locally-stored field data for a layer in a prismesh.
struct prismesh_layer_data_t
{
  /// Data storage for the layer. Use DECLARE_PRISMESH_*_ARRAY to provide 
  /// multidimensional array access to this data.
  void* data;
  
  /// The number of data indices within the x-y plane.
  size_t xy_size;

  /// The number of data indices along the z axis.
  size_t z_size;

  /// The centering of the data.
  prismesh_centering_t centering;

  /// The number of components for a datum in the field.
  size_t num_components;
};
typedef struct prismesh_layer_data_t prismesh_layer_data_t;

/// Constructs a new prismesh field with the given number of components
/// on the given mesh.
/// \memberof prismesh_field
prismesh_field_t* prismesh_field_new(prismesh_t* mesh,
                                     prismesh_centering_t centering,
                                     size_t num_components);

/// Destroys the given prismesh field.
/// \memberof prismesh_field
void prismesh_field_free(prismesh_field_t* field);

/// Returns the centering for the field.
/// \memberof prismesh_field
prismesh_centering_t prismesh_field_centering(prismesh_field_t* field);

/// Returns the number of components in the field.
/// \memberof prismesh_field
size_t prismesh_field_num_components(prismesh_field_t* field);

/// Returns the number of (locally stored) layers in the field.
/// \memberof prismesh_field
size_t prismesh_field_num_layers(prismesh_field_t* field);

/// Returns an internal pointer to the field's underlying prismesh.
/// \memberof prismesh_field
prismesh_t* prismesh_field_mesh(prismesh_field_t* field);

/// Traverses locally-stored field data layer by layer.
/// \param [in,out] pos Controls the traversal. Set to 0 to reset.
/// \param [out] layer Stores the layer found in the traversal.
/// \param [out] data Stores the layer data associated with \ref layer.
/// \returns true if a locally-stored layer is found, false if not.
/// \memberof prismesh_field
bool prismesh_field_next_layer(prismesh_field_t* field, int* pos, 
                               prismesh_layer_t** layer, 
                               prismesh_layer_data_t** data);

///@}

///@{
// These macros generate multidimensional arrays that can access the given
// layer data using C99 variable-length arrays.

/// \def DECLARE_PRISMESH_CELL_ARRAY
/// Allows access to prismesh cell data. Cell arrays are indexed the following 
/// way:
/// array[i][k][c] where i is the column index, k is the z index,
/// and c is the component.
/// * i runs from 0 to data->num_xy_data-1 for interior cells, with ghost 
///   cells beginning at data->num_xy_data and proceeding hence.
/// * k runs from 1 to data->num_vertical_cells for interior cells, with ghost
///   values at 0 and data->num_vertical_cells+1.
/// * c runs from 0 to data->num_components-1.
#define DECLARE_PRISMESH_CELL_ARRAY(array, data) \
ASSERT(data->centering == PRISMESH_CELL); \
DECLARE_3D_ARRAY(real_t, array, data->data, data->num_xy_data, data->num_vertical_cells, data->num_components)

/// \def DECLARE_PRISMESH_XYFACE_ARRAY
/// Allows access to prismesh data for faces with normal vectors in the x-y 
/// plane. XY-face arrays are indexed the following way:
/// array[i][k][c] with (i, k) identifying the ith xy-face at the kth 
/// vertical position, and c as the component. 
/// * i runs from 0 to data->num_xy_data.
/// * k runs from 0 to data->num_vertical_cells-1.
/// * c runs from 0 to data->num_components-1.
/// Use layer->cell_xyfaces[m][n] to retrieve the index for the nth xy-face 
/// of the mth cell.
#define DECLARE_PRISMESH_XYFACE_ARRAY(array, data) \
ASSERT(data->centering == PRISMESH_XYFACE); \
DECLARE_3D_ARRAY(real_t, array, data->data, data->num_xy_data, data->num_vertical_cells, data->num_components)

/// \def DECLARE_PRISMESH_ZFACE_ARRAY
/// Allows access to prismesh data for faces with normal vectors in the 
/// +/- z direction. Z-face arrays are indexed the following way:
/// array[i][k][c] identifies the "bottom" z-face for the ith cell at 
/// vertical level k, with c as the component.
/// * i runs from 0 to data->num_xy_data-1.
/// * k runs from 0 to data->num_vertical_cells.
/// * c runs from 0 to data->num_components-1.
#define DECLARE_PRISMESH_ZFACE_ARRAY(array, data) \
ASSERT(data->centering == PRISMESH_ZFACE); \
DECLARE_3D_ARRAY(real_t, array, data->data, data->num_xy_data, data->num_vertical_cells+1, data->num_components)

/// \def DECLARE_PRISMESH_XYEDGE_ARRAY
/// Allows access to prismesh edge data for edges that lie in the x-y plane. 
/// XY-edge arrays are indexed the following way:
/// array[i][k][c] identifies the ith xy-edge in the set of xy-edges aligned 
/// with the bottom of cells in the kth vertical position, with c as the 
/// component.
/// * i runs from 0 to data->num_xy_data-1.
/// * k runs from 0 to data->num_vertical_cells.
/// * c runs from 0 to data->num_components-1.
/// Use layer->zface_xyedges[m][n] to retrieve the index for the nth xy-edge 
/// of the mth z-face.
#define DECLARE_PRISMESH_XYEDGE_ARRAY(array, data) \
ASSERT(data->centering == PRISMESH_XYEDGE); \
DECLARE_3D_ARRAY(real_t, array, data->data, data->num_xy_data, data->num_vertical_cells+1, data->num_components)

/// \def DECLARE_PRISMESH_ZEDGE_ARRAY
/// Allows access to prismesh edge data for edges that align with the z axis.
/// Z-edge arrays are indexed the following way:
/// array[i][k][c] identifies the ith z-edge for xy-faces in the kth vertical
/// position, with c as the component.
/// * i runs from 0 to data->num_xy_data-1.
/// * k runs from 0 to data->num_vertical_cells-1.
/// * c runs from 0 to data->num_components-1.
/// Use layer->xyface_zedges[m][n] to retrieve the index for the nth z-edge 
/// of the mth xy-face.
#define DECLARE_PRISMESH_ZEDGE_ARRAY(array, data) \
ASSERT(data->centering == PRISMESH_ZEDGE); \
DECLARE_3D_ARRAY(real_t, array, data->data, data->num_xy_data, data->num_vertical_cells, data->num_components)

/// \def DECLARE_PRISMESH_NODE_ARRAY
/// Allows access to prismesh node data. Node arrays are indexed the 
/// following way:
/// array[i][k][c] identifies the ith node of the bottom z-face for the cell
/// in the kth vertical position, with c as the component.
/// * i runs from 0 to data->num_xy_data-1.
/// * k runs from 0 to data->num_vertical_cells.
/// * c runs from 0 to data->num_components-1.
/// Use layer->xyface_nodes[m][n] to retrieve the index for the nth node
/// of the mth xy-face, and layer->zface_nodes[m][n] to retrieve the index
/// for the nth node of the mth z-face.
#define DECLARE_PRISMESH_NODE_ARRAY(array, data) \
ASSERT(data->centering == PRISMESH_NODE); \
DECLARE_3D_ARRAY(real_t, array, data->data, data->num_xy_data, data->num_vertical_cells+1, data->num_components)

///@}

#endif

