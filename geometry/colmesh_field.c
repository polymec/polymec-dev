// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/enumerable.h"
#include "core/exchanger.h"
#include "core/timer.h"
#include "core/unordered_map.h"
#include "geometry/colmesh_field.h"

static colmesh_chunk_data_t* colmesh_chunk_data_with_buffer(colmesh_chunk_t* chunk,
                                                            colmesh_centering_t centering,
                                                            int num_components, 
                                                            void* buffer)
{
  colmesh_chunk_data_t* data = polymec_malloc(sizeof(colmesh_chunk_data_t));
  data->chunk = chunk;
  data->centering = centering;
  data->num_components = num_components;

  // Figure out the dimensions of our data.
  if (centering == COLMESH_CELL)
  {
    data->xy_size = chunk->num_columns + chunk->num_ghost_columns;
    data->z_size = chunk->num_z_cells + 2;
  }
  else if (centering == COLMESH_XYFACE)
  {
    data->xy_size = chunk->num_xy_faces;
    data->z_size = chunk->num_z_cells;
  }
  else if (centering == COLMESH_ZFACE)
  {
    data->xy_size = chunk->num_columns;
    data->z_size = chunk->num_z_cells + 1;
  }
  else if (centering == COLMESH_XYEDGE)
  {
    data->xy_size = chunk->num_xy_faces;
    data->z_size = chunk->num_z_cells + 1;
  }
  else if (centering == COLMESH_ZEDGE)
  {
    data->xy_size = chunk->num_xy_nodes;
    data->z_size = chunk->num_z_cells;
  }
  else // (centering == COLMESH_NODE)
  {
    data->xy_size = chunk->num_xy_nodes;
    data->z_size = chunk->num_z_cells + 1;
  }
  data->data = buffer;
  return data;
}

static void colmesh_chunk_data_free(colmesh_chunk_data_t* data)
{
  polymec_free(data);
}

void colmesh_chunk_data_copy(colmesh_chunk_data_t* data, 
                             colmesh_chunk_data_t* dest)
{
  ASSERT(dest->chunk->num_columns == data->chunk->num_columns);
  ASSERT(dest->chunk->num_z_cells == data->chunk->num_z_cells);
  ASSERT(dest->centering == data->centering);
  ASSERT(dest->num_components == data->num_components);
  ASSERT(dest->xy_size == data->xy_size);
  ASSERT(dest->z_size == data->z_size);
  memcpy(dest->data, data->data, colmesh_chunk_data_size(data));
}

DEFINE_UNORDERED_MAP(chunk_data_map, int, colmesh_chunk_data_t*, int_hash, int_equals)

struct colmesh_field_t 
{
  colmesh_t* mesh;
  colmesh_centering_t centering;
  int num_components;

  // Chunks and data storage.
  size_t num_xy_chunks, num_z_chunks;
  chunk_data_map_t* chunks;
  void* buffer;
  size_t bytes;
  bool owns_buffer;

  // Exchangers.
  exchanger_t *ex;
  int ex_token;

  // Metadata.
  field_metadata_t* md;
};

static inline int chunk_index(colmesh_field_t* field, int xy_index, int z_index)
{
  return (int)(field->num_z_chunks * xy_index + z_index);
}

extern exchanger_t* colmesh_exchanger(colmesh_t* mesh, 
                                      colmesh_centering_t centering);

colmesh_field_t* colmesh_field_with_buffer(colmesh_t* mesh,
                                           colmesh_centering_t centering,
                                           int num_components,
                                           void* buffer,
                                           bool assume_control)
{
  START_FUNCTION_TIMER();
  ASSERT(num_components > 0);
  colmesh_field_t* field = polymec_malloc(sizeof(colmesh_field_t));
  field->mesh = mesh;
  field->centering = centering;
  field->num_components = num_components;

  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  colmesh_get_chunk_info(mesh, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  field->num_xy_chunks = num_xy_chunks;
  field->num_z_chunks = num_z_chunks;

  field->chunks = chunk_data_map_new();
  field->buffer = NULL;
  field->bytes = 0;
  field->owns_buffer = false;

  // Set up exchangers and tokens.
  field->ex = colmesh_exchanger(mesh, centering);
  if (field->ex != NULL) // FIXME: Remove this test when all exchangers are complete!
    retain_ref(field->ex);
  field->ex_token = -1;

  // Now populate the chunks (with NULL buffers).
  int pos = 0, xy_index, z_index;
  colmesh_chunk_t* chunk;
  while (colmesh_next_chunk(mesh, &pos, &xy_index, &z_index, &chunk))
  {
    int ch_index = chunk_index(field, xy_index, z_index);
    colmesh_chunk_data_t* chunk_data = colmesh_chunk_data_with_buffer(chunk, centering, num_components, NULL);
    chunk_data_map_insert_with_v_dtor(field->chunks, ch_index, chunk_data, colmesh_chunk_data_free);
    field->bytes += colmesh_chunk_data_size(chunk_data);
  }

  // Use the given buffer.
  colmesh_field_set_buffer(field, buffer, assume_control);

  // Set up metadata.
  field->md = field_metadata_new(num_components);

  STOP_FUNCTION_TIMER();
  return field;
}

colmesh_field_t* colmesh_field_new(colmesh_t* mesh,
                                   colmesh_centering_t centering,
                                   int num_components)
{
  ASSERT(num_components > 0);
  START_FUNCTION_TIMER();
  colmesh_field_t* field = colmesh_field_with_buffer(mesh, centering, num_components, NULL, false);
  void* buffer = polymec_calloc(field->bytes, 1);
  colmesh_field_set_buffer(field, buffer, true);
  STOP_FUNCTION_TIMER();
  return field;
}

void colmesh_field_free(colmesh_field_t* field)
{
  chunk_data_map_free(field->chunks);
  if ((field->buffer != NULL) && field->owns_buffer)
    polymec_free(field->buffer);
  if (field->ex != NULL) 
    release_ref(field->ex);
  release_ref(field->md);
  polymec_free(field);
}

colmesh_centering_t colmesh_field_centering(colmesh_field_t* field)
{
  return field->centering;
}

int colmesh_field_num_components(colmesh_field_t* field)
{
  return field->num_components;
}

size_t colmesh_field_num_chunks(colmesh_field_t* field)
{
  return field->chunks->size;
}

colmesh_t* colmesh_field_mesh(colmesh_field_t* field)
{
  return field->mesh;
}

void* colmesh_field_buffer(colmesh_field_t* field)
{
  return field->buffer;
}

void colmesh_field_set_buffer(colmesh_field_t* field, 
                              void* buffer, 
                              bool assume_control)
{
  START_FUNCTION_TIMER();
  if ((field->buffer != NULL) && field->owns_buffer)
    polymec_free(field->buffer);
  field->buffer = buffer;
  field->owns_buffer = assume_control;

  // Point the chunks at the buffer.
  int pos = 0, xy_index, z_index;
  colmesh_chunk_data_t* chunk_data;
  size_t offset = 0;
  while (colmesh_field_next_chunk(field, &pos, &xy_index, &z_index, &chunk_data))
  {
    chunk_data->data = &(((real_t*)buffer)[offset]);
    offset += colmesh_chunk_data_size(chunk_data) / sizeof(real_t);
  }
  ASSERT(offset == (field->bytes / sizeof(real_t)));
  STOP_FUNCTION_TIMER();
}

void colmesh_field_copy(colmesh_field_t* field,
                        colmesh_field_t* dest)
{
  ASSERT(dest->mesh == field->mesh);
  ASSERT(dest->centering == field->centering);
  ASSERT(dest->num_components == field->num_components);
  ASSERT(dest->num_xy_chunks == field->num_xy_chunks);
  ASSERT(dest->num_z_chunks == field->num_z_chunks);
  START_FUNCTION_TIMER();
  int pos = 0, ch_index;
  colmesh_chunk_data_t* src_data;
  while (chunk_data_map_next(field->chunks, &pos, &ch_index, &src_data))
  {
    colmesh_chunk_data_t** dest_data_p = chunk_data_map_get(dest->chunks, ch_index);
    ASSERT(dest_data_p != NULL);
    colmesh_chunk_data_t* dest_data = *dest_data_p;
    ASSERT(dest_data->chunk == src_data->chunk);
    colmesh_chunk_data_copy(src_data, dest_data);
  }

  // Copy metadata.
  release_ref(dest->md);
  dest->md = field_metadata_clone(field->md);

  STOP_FUNCTION_TIMER();
}

field_metadata_t* colmesh_field_metadata(colmesh_field_t* field)
{
  return field->md;
}

colmesh_chunk_data_t* colmesh_field_chunk_data(colmesh_field_t* field, 
                                               int xy_index,
                                               int z_index)
{
  int index = chunk_index(field, xy_index, z_index);
  colmesh_chunk_data_t** data_p = chunk_data_map_get(field->chunks, index);
  if (data_p != NULL)
    return *data_p;
  else
    return NULL;
}

bool colmesh_field_next_chunk(colmesh_field_t* field, int* pos, 
                              int* xy_index, int* z_index,
                              colmesh_chunk_data_t** chunk_data)
{
  colmesh_chunk_t* chunk;
  bool result = colmesh_next_chunk(field->mesh, pos, xy_index, z_index, &chunk);
  if (result)
  {
    *chunk_data = colmesh_field_chunk_data(field, *xy_index, *z_index);
    ASSERT((*chunk_data)->chunk == chunk);
  }
  return result;
}

void colmesh_field_exchange(colmesh_field_t* field)
{
  colmesh_field_start_exchange(field);
  colmesh_field_finish_exchange(field);
}

void colmesh_field_start_exchange(colmesh_field_t* field)
{
  ASSERT(!colmesh_field_is_exchanging(field));
  START_FUNCTION_TIMER();

  // Do we have an exchanger yet?
  if (field->ex == NULL)
  {
    field->ex = colmesh_exchanger(field->mesh, field->centering);
    retain_ref(field->ex);
  }

  // Start the xy exchange.
  int stride = (int)field->num_components;
  field->ex_token = exchanger_start_exchange(field->ex, field->buffer, stride, 0, MPI_REAL_T);
  STOP_FUNCTION_TIMER();
}

void colmesh_field_finish_exchange(colmesh_field_t* field)
{
  ASSERT(colmesh_field_is_exchanging(field));
  START_FUNCTION_TIMER();
  if (field->ex_token != -1)
    exchanger_finish_exchange(field->ex, field->ex_token);
  field->ex_token = -1;
  STOP_FUNCTION_TIMER();
}

bool colmesh_field_is_exchanging(colmesh_field_t* field)
{
  return (field->ex_token != -1);
}

void colmesh_field_set_exchanger(colmesh_field_t* field, exchanger_t* ex)
{
  if (field->ex != NULL)
    release_ref(field->ex);
  if (ex != NULL)
    retain_ref(ex);
  field->ex = ex;
}


real_enumerable_generator_t* colmesh_field_enumerate(colmesh_field_t* field)
{
  size_t num_values = 0;
  int pos = 0, xy, z;
  colmesh_chunk_data_t* chunk_data;
  while (colmesh_field_next_chunk(field, &pos, &xy, &z, &chunk_data))
    num_values += colmesh_chunk_data_size(chunk_data) / sizeof(real_t);

  size_t offset = 0;
  real_t* array = polymec_malloc(sizeof(real_t) * num_values);
  pos = 0;
  while (colmesh_field_next_chunk(field, &pos, &xy, &z, &chunk_data))
  {
    size_t data_size = colmesh_chunk_data_size(chunk_data);
    memcpy(&array[offset], chunk_data->data, data_size);
    offset += data_size / sizeof(real_t);
  }
  return real_enumerable_generator_from_array(array, num_values, true);
}

real_enumerable_generator_t* colmesh_chunk_data_enumerate(colmesh_chunk_data_t* chunk_data)
{
  size_t num_values = colmesh_chunk_data_size(chunk_data) / sizeof(real_t);
  return real_enumerable_generator_from_array((real_t*)chunk_data->data, num_values, NULL);
}

