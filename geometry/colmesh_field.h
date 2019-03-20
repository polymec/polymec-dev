// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_COLMESH_FIELD_H
#define POLYMEC_COLMESH_FIELD_H

#include "core/declare_nd_array.h"
#include "geometry/colmesh.h"
#include "geometry/field_metadata.h"

/// \addtogroup geometry geometry
///@{

/// \class colmesh_field
/// This type represents field values defined on a colmesh.
typedef struct colmesh_field_t colmesh_field_t;

/// \struct colmesh_chunk_data
/// This type holds locally-stored field data for a chunk in a colmesh.
struct colmesh_chunk_data_t
{
  /// The chunk for this chunk data. Provided for convenience.
  colmesh_chunk_t* chunk;

  /// Data storage for the chunk. Use DECLARE_COLMESH_*_ARRAY to provide
  /// multidimensional array access to this data.
  void* data;

  /// The number of data indices within the x-y plane. This corresponds to
  /// a quantity in the chunk, depending on the centering.
  size_t xy_size;

  /// The number of data indices along the z axis.
  size_t z_size;

  /// The centering of the data.
  colmesh_centering_t centering;

  /// The number of components for a datum in the field.
  int num_components;
};
typedef struct colmesh_chunk_data_t colmesh_chunk_data_t;

/// Returns the number of bytes occupied by the chunk data.
/// \memberof colmesh_chunk_data
static inline size_t colmesh_chunk_data_size(colmesh_chunk_data_t* data)
{
  return sizeof(real_t) * data->xy_size * data->z_size * data->num_components;
}

/// Copies the data from this chunk to a destination chunk.
/// \memberof colmesh_chunk_data
void colmesh_chunk_data_copy(colmesh_chunk_data_t* data,
                             colmesh_chunk_data_t* dest);

/// Constructs a new colmesh field with the given number of components
/// on the given mesh.
/// \param [in] mesh The colmesh on which the field is created. Must be finalized.
/// \param [in] centering The centering of which the new field.
/// \param [in] num_components The number of components in the new field.
/// \memberof colmesh_field
colmesh_field_t* colmesh_field_new(colmesh_t* mesh,
                                   colmesh_centering_t centering,
                                   int num_components);

/// Constructs a new colmesh field that uses the given buffer for data storage.
/// \param [in] mesh The colmesh on which the field is created. Must be finalized.
/// \param [in] centering The centering of which the new field.
/// \param [in] num_components The number of components in the new field.
/// \param [in] buffer The buffer to use for storing field data. Must be ample for storing all the
///                    field's data.
/// \param [in] assume_control If true, the field assumes ownership of the
///                            buffer. Otherwise it does not.
/// \memberof colmesh_field
colmesh_field_t* colmesh_field_with_buffer(colmesh_t* mesh,
                                           colmesh_centering_t centering,
                                           int num_components,
                                           void* buffer,
                                           bool assume_control);

/// Destroys the given colmesh field.
/// \memberof colmesh_field
void colmesh_field_free(colmesh_field_t* field);

/// Returns the buffer used to store the field's data.
/// \memberof colmesh_field
void* colmesh_field_buffer(colmesh_field_t* field);

/// Sets the buffer used by the field to store its data.
/// \param [in] buffer The buffer used by the field. Must be ample for storing the field's data.
/// \param [in] assume_control If true, the field assumes ownership of the buffer. Otherwise it does not.
/// \memberof colmesh_field
void colmesh_field_set_buffer(colmesh_field_t* field,
                              void* buffer,
                              bool assume_control);

/// Copies the data in this field to a destination field.
/// \memberof colmesh_field
void colmesh_field_copy(colmesh_field_t* field, colmesh_field_t* dest);

/// Returns the metadata associated with this field. Every field has a
/// metadata object that is empty until its properties are specified.
/// \memberof unimesh_field
field_metadata_t* colmesh_field_metadata(colmesh_field_t* field);

/// Returns the centering for the field.
/// \memberof colmesh_field
colmesh_centering_t colmesh_field_centering(colmesh_field_t* field);

/// Returns the number of components in the field.
/// \memberof colmesh_field
int colmesh_field_num_components(colmesh_field_t* field);

/// Returns the number of (locally stored) chunks in the field.
/// \memberof colmesh_field
size_t colmesh_field_num_chunks(colmesh_field_t* field);

/// Returns an internal pointer to the field's underlying colmesh.
/// \memberof colmesh_field
colmesh_t* colmesh_field_mesh(colmesh_field_t* field);

/// Returns the chunk data for the chunk in this field identified by the
/// given xy and z indices.
/// \memberof colmesh_field
colmesh_chunk_data_t* colmesh_field_chunk_data(colmesh_field_t* field,
                                               int xy_index,
                                               int z_index);

/// Traverses locally-stored field data chunk by chunk.
/// \param pos [in,out] Controls the traversal. Set to 0 to reset.
/// \param xy_index [out] Stores the xy index of the next chunk.
/// \param z_index [out] Stores the z index of the next chunk.
/// \param chunk_data [out] Stores the chunk data for the chunk. Note that
///                         the chunk itself is available via chunk_data->chunk.
/// \returns true if a locally-stored chunk is found, false if not.
/// \memberof colmesh_field
bool colmesh_field_next_chunk(colmesh_field_t* field, int* pos,
                              int* xy_index, int* z_index,
                              colmesh_chunk_data_t** chunk_data);

/// Synchronously exchanges boundary data in the chunks within this
/// field with that of adjoining chunks. For cell-centered data, this
/// means filling ghost cells. For face-, node-, and edge-centered data, it
/// means overwriting values on the boundary of each chunk with data from
/// other chunks.
/// \memberof colmesh_field
void colmesh_field_exchange(colmesh_field_t* field);

/// Begins an asynchronous exchange of boundary data for the chunks in this
/// field.
/// \memberof colmesh_field
void colmesh_field_start_exchange(colmesh_field_t* field);

/// Finishes an asynchronous exchange initiated with
/// \ref colmesh_field_start_exchange.
/// \memberof colmesh_field
void colmesh_field_finish_exchange(colmesh_field_t* field);

/// Returns `true` if this field is in the middle of an asynchronous exchange,
/// `false` if not.
/// \memberof colmesh_field
bool colmesh_field_is_exchanging(colmesh_field_t* field);

/// Sets the exchanger used by the field for exchanges.
/// Use this method instead of directly assigning a new exchanger to the field.
/// \param [in] ex The exchanger to be used by this field.
/// \memberof colmesh_field
void colmesh_field_set_exchanger(colmesh_field_t* field, exchanger_t* ex);

typedef struct real_enumerable_generator_t real_enumerable_generator_t;

/// Enumerates values in the given colmesh field.
/// \memberof colmesh_field
real_enumerable_generator_t* colmesh_field_enumerate(colmesh_field_t* field);

/// Enumerates values in the given colmesh chunk data set.
/// \memberof colmesh_chunk_data
real_enumerable_generator_t* colmesh_chunk_data_enumerate(colmesh_chunk_data_t* chunk_data);

///@}

///@{
// These macros generate multidimensional arrays that can access the given
// chunk data using C99 variable-length arrays.

/// \def DECLARE_COLMESH_CELL_ARRAY(array, chunk_data)
/// Allows access to colmesh cell data. Cell arrays are indexed the following
/// way:
/// `array[xy][z][c]` where `i` is the xy index, `z` is the z index,
/// and `c` is the component.
/// * `xy` runs from 0 to `chunk_data->chunk->num_columns-1` for interior cells,
///   with ghost cells beginning at `chunk->num_columns` and proceeding hence.
/// * `z` runs from 1 to `chunk_data->chunk->num_z_cells` for interior cells, with
///   ghost values at 0 and `chunk_data->chunk->num_z_cells+1`.
/// * `c` runs from 0 to chunk->num_components-1.
#define DECLARE_COLMESH_CELL_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == COLMESH_CELL); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_COLMESH_XYFACE_ARRAY(array, chunk_data)
/// Allows access to colmesh data for faces with normal vectors in the xy
/// plane. XY-face arrays are indexed the following way:
/// `array[xy][z][c]` with (xy, z) identifying a face at the given xy and z position,
/// and `c` as the component.
/// * `xy` runs from 0 to `chunk_data->chunk->num_xy_faces-1`.
/// * `z` runs from 0 to `chunk_data->chunk->num_z_cells-1`.
/// * `c` runs from 0 to `chunk_data->num_components-1`.
#define DECLARE_COLMESH_XYFACE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == COLMESH_XYFACE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_COLMESH_ZFACE_ARRAY(array, chunk_data)
/// Allows access to colmesh data for faces with normal vectors in the
/// +/- z direction. Z-face arrays are indexed the following way:
/// `array[xy][z][c]` identifies the "bottom" z-face for the cell at (xy, z),
/// with c as the component.
/// * `xy` runs from 0 to `chunk_data->chunk->num_columns-1`.
/// * `z` runs from 0 to `chunk_data->chunk->num_z_cells`.
/// * `c` runs from 0 to `chunk_data->num_components-1`.
#define DECLARE_COLMESH_ZFACE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == COLMESH_ZFACE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_COLMESH_XYEDGE_ARRAY(array, chunk_data)
/// Allows access to colmesh edge data for edges that lie in the xy plane.
/// XY-edge arrays are indexed the following way:
/// `array[xy][k][c]` identifies the xy-edge at the given xy position, running
/// along the bottom of cells in the given z position, with `c` as the
/// component.
/// * `xy` runs from 0 to `chunk_data->chunk->num_xy_edges-1`.
/// * `z` runs from 0 to `chunk_data->chunk->num_z_cells`.
/// * `c` runs from 0 to `chunk_data->num_components-1`.
#define DECLARE_COLMESH_XYEDGE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == COLMESH_XYEDGE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_COLMESH_ZEDGE_ARRAY(array, chunk_data)
/// Allows access to colmesh edge data for edges that align with the z axis.
/// z-edge arrays are indexed the following way:
/// `array[xy][z][c]` identifies the z-edge at the given xy, running vertically along
/// cells in the given z position, with c as the component.
/// * `xy` runs from 0 to `chunk_data->chunk->num_xy_nodes-1`.
/// * `k` runs from 0 to `chunk_data->chunk->num_z_cells-1`.
/// * `c` runs from 0 to `chunk_data->num_components-1`.
#define DECLARE_COLMESH_ZEDGE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == COLMESH_ZEDGE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

/// \def DECLARE_COLMESH_NODE_ARRAY(array, chunk_data)
/// Allows access to colmesh node data. Node arrays are indexed the
/// following way:
/// `array[xy][k][c]` identifies the node at the given xy position along the bottom
/// z-face for the cell at the given z position, with c as the component.
/// * `xy` runs from 0 to `chunk_data->chunk->num_xy_nodes-1`.
/// * `z` runs from 0 to `chunk_data->chunk->num_z_cells`.
/// * `c` runs from 0 to `chunk->num_components-1`.
#define DECLARE_COLMESH_NODE_ARRAY(array, chunk_data) \
ASSERT(chunk_data->centering == COLMESH_NODE); \
DECLARE_3D_ARRAY(real_t, array, chunk_data->data, chunk_data->xy_size, chunk_data->z_size, chunk_data->num_components)

///@}

#endif

