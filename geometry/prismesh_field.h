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

/// \struct prismesh_chunk_data
/// This type holds locally-stored field data for a chunk in a prismesh.
struct prismesh_chunk_data_t
{
  /// Data storage for the chunk. Use DECLARE_PRISMESH_*_ARRAY to provide 
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
typedef struct prismesh_chunk_data_t prismesh_chunk_data_t;

/// Returns the number of bytes occupied by the chunk data.
/// \memberof prismesh_chunk_data
static inline size_t prismesh_chunk_data_size(prismesh_chunk_data_t* data)
{
  return sizeof(real_t) * data->xy_size * data->z_size * data->num_components;
}

/// Copies the data from this chunk to a destination chunk.
/// \memberof prismesh_chunk_data
void prismesh_chunk_data_copy(prismesh_chunk_data_t* data, 
                              prismesh_chunk_data_t* dest);

/// Constructs a new prismesh field with the given number of components
/// on the given mesh.
/// \memberof prismesh_field
prismesh_field_t* prismesh_field_new(prismesh_t* mesh,
                                     prismesh_centering_t centering,
                                     size_t num_components);

/// Destroys the given prismesh field.
/// \memberof prismesh_field
void prismesh_field_free(prismesh_field_t* field);

/// Copies the data in this field to a destination field.
/// \memberof prismesh_field
void prismesh_field_copy(prismesh_field_t* field,
                         prismesh_field_t* dest);

/// Returns the centering for the field.
/// \memberof prismesh_field
prismesh_centering_t prismesh_field_centering(prismesh_field_t* field);

/// Returns the number of components in the field.
/// \memberof prismesh_field
size_t prismesh_field_num_components(prismesh_field_t* field);

/// Returns the number of (locally stored) chunks in the field.
/// \memberof prismesh_field
size_t prismesh_field_num_chunks(prismesh_field_t* field);

/// Returns an internal pointer to the field's underlying prismesh.
/// \memberof prismesh_field
prismesh_t* prismesh_field_mesh(prismesh_field_t* field);

/// Returns the chunk data for the chunk in this field identified by the 
/// given xy and z indices.
/// \memberof prismesh_field
prismesh_chunk_data_t* prismesh_field_chunk_data(prismesh_field_t* field, 
                                                 int xy_index,
                                                 int z_index);

/// Traverses locally-stored field data chunk by chunk.
/// \param pos [in,out] Controls the traversal. Set to 0 to reset.
/// \param xy_index [out] Stores the xy index of the next chunk.
/// \param z_index [out] Stores the z index of the next chunk.
/// \param chunk [out] Stores the chunk found in the traversal.
/// \param data [out] Stores the chunk data associated with \ref chunk.
/// \returns true if a locally-stored chunk is found, false if not.
/// \memberof prismesh_field
bool prismesh_field_next_chunk(prismesh_field_t* field, int* pos, 
                               int* xy_index, int* z_index,
                               prismesh_chunk_t** chunk, 
                               prismesh_chunk_data_t** data);

///@}

///@{
// These macros generate multidimensional arrays that can access the given
// chunk data using C99 variable-length arrays.

/// \def DECLARE_PRISMESH_CELL_ARRAY
/// Allows access to prismesh cell data. Cell arrays are indexed the following 
/// way:
/// array[i][k][c] where i is the column index, k is the z index,
/// and c is the component.
/// * i runs from 0 to chunk->xy_size-1 for interior cells, with ghost 
///   cells beginning at chunk->xy_size and proceeding hence.
/// * k runs from 1 to chunk->z_size for interior cells, with ghost
///   values at 0 and chunk->z_size+1.
/// * c runs from 0 to chunk->num_components-1.
#define DECLARE_PRISMESH_CELL_ARRAY(array, chunk) \
ASSERT(chunk->centering == PRISMESH_CELL); \
DECLARE_3D_ARRAY(real_t, array, chunk->data, chunk->xy_size, chunk->z_size+2, chunk->num_components)

/// \def DECLARE_PRISMESH_XYFACE_ARRAY
/// Allows access to prismesh data for faces with normal vectors in the x-y 
/// plane. XY-face arrays are indexed the following way:
/// array[i][k][c] with (i, k) identifying the ith xy-face at the kth 
/// vertical position, and c as the component. 
/// * i runs from 0 to chunk->xy_size-1.
/// * k runs from 0 to chunk->z_size-1.
/// * c runs from 0 to chunk->num_components-1.
/// Use chunk->cell_xy_faces[m][n] to retrieve the index for the nth xy-face 
/// of the mth cell.
#define DECLARE_PRISMESH_XYFACE_ARRAY(array, chunk) \
ASSERT(chunk->centering == PRISMESH_XYFACE); \
DECLARE_3D_ARRAY(real_t, array, chunk->data, chunk->xy_size, chunk->z_size, chunk->num_components)

/// \def DECLARE_PRISMESH_ZFACE_ARRAY
/// Allows access to prismesh data for faces with normal vectors in the 
/// +/- z direction. Z-face arrays are indexed the following way:
/// array[i][k][c] identifies the "bottom" z-face for the ith cell at 
/// vertical level k, with c as the component.
/// * i runs from 0 to chunk->xy_size-1.
/// * k runs from 0 to chunk->z_size.
/// * c runs from 0 to chunk->num_components-1.
#define DECLARE_PRISMESH_ZFACE_ARRAY(array, chunk) \
ASSERT(chunk->centering == PRISMESH_ZFACE); \
DECLARE_3D_ARRAY(real_t, array, chunk->data, chunk->xy_size, chunk->z_size+1, chunk->num_components)

/// \def DECLARE_PRISMESH_XYEDGE_ARRAY
/// Allows access to prismesh edge data for edges that lie in the x-y plane. 
/// XY-edge arrays are indexed the following way:
/// array[i][k][c] identifies the ith xy-edge in the set of xy-edges aligned 
/// with the bottom of cells in the kth vertical position, with c as the 
/// component.
/// * i runs from 0 to chunk->xy_size-1.
/// * k runs from 0 to chunk->z_size.
/// * c runs from 0 to chunk->num_components-1.
/// Use chunk->zface_xyedges[m][n] to retrieve the index for the nth xy-edge 
/// of the mth z-face.
#define DECLARE_PRISMESH_XYEDGE_ARRAY(array, chunk) \
ASSERT(chunk->centering == PRISMESH_XYEDGE); \
DECLARE_3D_ARRAY(real_t, array, chunk->data, chunk->xy_size, chunk->z_size+1, chunk->num_components)

/// \def DECLARE_PRISMESH_ZEDGE_ARRAY
/// Allows access to prismesh edge data for edges that align with the z axis.
/// Z-edge arrays are indexed the following way:
/// array[i][k][c] identifies the ith z-edge for xy-faces in the kth vertical
/// position, with c as the component.
/// * i runs from 0 to chunk->xy_size-1.
/// * k runs from 0 to chunk->z_size-1.
/// * c runs from 0 to chunk->num_components-1.
/// Z-edges are indexed the same as nodes, so use chunk->xy_face_nodes[m][n] to 
/// retrieve the index for the nth z-edge of the mth xy-face.
#define DECLARE_PRISMESH_ZEDGE_ARRAY(array, chunk) \
ASSERT(chunk->centering == PRISMESH_ZEDGE); \
DECLARE_3D_ARRAY(real_t, array, chunk->data, chunk->xy_size, chunk->z_size, chunk->num_components)

/// \def DECLARE_PRISMESH_NODE_ARRAY
/// Allows access to prismesh node data. Node arrays are indexed the 
/// following way:
/// array[i][k][c] identifies the ith node of the bottom z-face for the cell
/// in the kth vertical position, with c as the component.
/// * i runs from 0 to chunk->xy_size-1.
/// * k runs from 0 to chunk->z_size.
/// * c runs from 0 to chunk->num_components-1.
/// Use chunk->xy_face_nodes[m][n] to retrieve the index for the nth node
/// of the mth xy-face, and chunk->zface_nodes[m][n] to retrieve the index
/// for the nth node of the mth z-face.
#define DECLARE_PRISMESH_NODE_ARRAY(array, chunk) \
ASSERT(chunk->centering == PRISMESH_NODE); \
DECLARE_3D_ARRAY(real_t, array, chunk->data, chunk->xy_size, chunk->z_size+1, chunk->num_components)

///@}

#endif

