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
  /// The chunk for this chunk data. Provided for convenience.
  prismesh_chunk_t* chunk;

  /// Data storage for the chunk. Use DECLARE_PRISMESH_*_ARRAY to provide 
  /// multidimensional array access to this data.
  void* data;
  
  /// The number of data indices within the x-y plane. This corresponds to 
  /// a quantity in the chunk, depending on the centering.
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
/// \param chunk_data [out] Stores the chunk data for the chunk. Note that 
///                         the chunk itself is available via chunk_data->chunk.
/// \returns true if a locally-stored chunk is found, false if not.
/// \memberof prismesh_field
bool prismesh_field_next_chunk(prismesh_field_t* field, int* pos, 
                               int* xy_index, int* z_index,
                               prismesh_chunk_data_t** chunk_data);

/// Compares all elements in the given component of the two given fields, 
/// returning true if the pairwise comparison of each component element in 
/// the two fields yields a "true" comparison, and false otherwise.
/// Calling this function on two fields with different centerings or with 
/// incompatible chunks is not allowed.
/// \memberof prismesh_field
bool prismesh_field_compare_all(prismesh_field_t* field,
                                prismesh_field_t* other_field,
                                int component,
                                bool (*comparator)(real_t val, real_t other_val));

/// Compares elements in the given component of the two given fields, 
/// returning true if the pairwise comparison of ANY component element in 
/// the two fields yields a "true" comparison, and false if none do so.
/// Calling this function on two fields with different centerings or with 
/// incompatible chunks is not allowed.
/// \memberof prismesh_field
bool prismesh_field_compare_any(prismesh_field_t* field,
                                prismesh_field_t* other_field,
                                int component,
                                bool (*comparator)(real_t val, real_t other_val));

/// Compares elements in the given component of the two given fields, 
/// returning true if NO pairwise comparison of ANY component element in 
/// the two fields yields a "true" comparison, and false if any do.
/// Calling this function on two fields with different centerings or with 
/// incompatible chunks is not allowed.
/// \memberof prismesh_field
bool prismesh_field_compare_none(prismesh_field_t* field,
                                 prismesh_field_t* other_field,
                                 int component,
                                 bool (*comparator)(real_t val, real_t other_val));

/// Compares all elements in the given component of the two given sets of chunk data, 
/// returning true if the pairwise comparison of each component element in 
/// the two datasets yields a "true" comparison, and false otherwise.
/// Calling this function on two datasets with different centerings or with 
/// incompatible chunks is not allowed.
/// \memberof prismesh_chunk_data
bool prismesh_chunk_data_compare_all(prismesh_chunk_data_t* chunk_data,
                                     prismesh_chunk_data_t* other_chunk_data,
                                     int component,
                                     bool (*comparator)(real_t val, real_t other_val));

/// Compares elements in the given component of the two given sets of chunk data, 
/// returning true if the pairwise comparison of ANY component element in 
/// the two datasets yields a "true" comparison, and false if none do so.
/// Calling this function on two datasets with different centerings or with 
/// incompatible chunks is not allowed.
/// \memberof prismesh_chunk_data
bool prismesh_chunk_data_compare_any(prismesh_chunk_data_t* chunk_data,
                                     prismesh_chunk_data_t* other_chunk_data,
                                     int component,
                                     bool (*comparator)(real_t val, real_t other_val));

/// Compares elements in the given component of the two given sets of chunk data, 
/// returning true if NO pairwise comparison of ANY component element in 
/// the two datasets yields a "true" comparison, and false if any do.
/// Calling this function on two datasets with different centerings or with 
/// incompatible chunks is not allowed.
/// \memberof prismesh_chunk_data
bool prismesh_chunk_data_compare_none(prismesh_chunk_data_t* chunk_data,
                                      prismesh_chunk_data_t* other_chunk_data,
                                      int component,
                                      bool (*comparator)(real_t val, real_t other_val));

///@}

///@{
// These macros generate multidimensional arrays that can access the given
// chunk data using C99 variable-length arrays.

/// \def DECLARE_PRISMESH_CELL_ARRAY
/// Allows access to prismesh cell data. Cell arrays are indexed the following 
/// way:
/// array[xy][z][c] where i is the xy index, z is the z index,
/// and c is the component.
/// * xy runs from 0 to chunk_data->chunk->num_columns-1 for interior cells, 
///   with ghost cells beginning at chunk->num_columns and proceeding hence.
/// * z runs from 1 to chunk_data->chunk->num_z_cells for interior cells, with 
///   ghost values at 0 and chunk->z_size+1.
/// * c runs from 0 to chunk->num_components-1.
#define DECLARE_PRISMESH_CELL_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == PRISMESH_CELL); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_PRISMESH_XYFACE_ARRAY
/// Allows access to prismesh data for faces with normal vectors in the xy 
/// plane. XY-face arrays are indexed the following way:
/// array[xy][z][c] with (xy, z) identifying a face at the given xy and z position,
/// and c as the component. 
/// * xy runs from 0 to chunk_data->chunk->num_xy_faces-1.
/// * z runs from 0 to chunk_data->chunk->num_z_cells-1.
/// * c runs from 0 to chunk_data->num_components-1.
#define DECLARE_PRISMESH_XYFACE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == PRISMESH_XYFACE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_PRISMESH_ZFACE_ARRAY
/// Allows access to prismesh data for faces with normal vectors in the 
/// +/- z direction. Z-face arrays are indexed the following way:
/// array[xy][z][c] identifies the "bottom" z-face for the cell at (xy, z),
/// with c as the component.
/// * xy runs from 0 to chunk_data->chunk->num_columns-1.
/// * z runs from 0 to chunk_data->chunk->num_z_cells.
/// * c runs from 0 to chunk_data->num_components-1.
#define DECLARE_PRISMESH_ZFACE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == PRISMESH_ZFACE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_PRISMESH_XYEDGE_ARRAY
/// Allows access to prismesh edge data for edges that lie in the xy plane. 
/// XY-edge arrays are indexed the following way:
/// array[xy][k][c] identifies the xy-edge at the given xy position, running
/// along the bottom of cells in the given z position, with c as the 
/// component.
/// * xy runs from 0 to chunk_data->chunk->num_xy_edges-1.
/// * z runs from 0 to chunk_data->chunk->num_z_cells.
/// * c runs from 0 to chunk_data->num_components-1.
#define DECLARE_PRISMESH_XYEDGE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == PRISMESH_XYEDGE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_PRISMESH_ZEDGE_ARRAY
/// Allows access to prismesh edge data for edges that align with the z axis.
/// z-edge arrays are indexed the following way:
/// array[xy][z][c] identifies the z-edge at the given xy, running vertically along
/// cells in the given z position, with c as the component.
/// * xy runs from 0 to chunk_data->chunk->num_xy_nodes-1.
/// * k runs from 0 to chunk_data->chunk->num_z_cells-1.
/// * c runs from 0 to chunk_data->num_components-1.
#define DECLARE_PRISMESH_ZEDGE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == PRISMESH_ZEDGE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_PRISMESH_NODE_ARRAY
/// Allows access to prismesh node data. Node arrays are indexed the 
/// following way:
/// array[xy][k][c] identifies the node at the given xy position along the bottom 
/// z-face for the cell at the given z position, with c as the component.
/// * xy runs from 0 to chunk_data->chunk->num_xy_nodes-1.
/// * z runs from 0 to chunk_data->chunk->num_z_cells.
/// * c runs from 0 to chunk->num_components-1.
#define DECLARE_PRISMESH_NODE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == PRISMESH_NODE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

///@}

#endif

