// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "geometry/prismesh_field.h"

static prismesh_chunk_data_t* prismesh_chunk_data_new(prismesh_chunk_t* chunk,
                                                      prismesh_centering_t centering,
                                                      size_t num_components)
{
  prismesh_chunk_data_t* data = polymec_malloc(sizeof(prismesh_chunk_data_t));
  data->chunk = chunk;
  data->centering = centering;
  data->num_components = num_components;

  // Figure out the dimensions of our data.
  if (centering == PRISMESH_CELL)
  {
    size_t num_ghost_xyvals = 0;
    data->xy_size = chunk->num_columns + num_ghost_xyvals;
    data->z_size = chunk->num_z_cells + 2;
  }
  else if (centering == PRISMESH_XYFACE)
  {
    data->xy_size = chunk->num_xy_faces;
    data->z_size = chunk->num_z_cells;
  }
  else if (centering == PRISMESH_ZFACE)
  {
    data->xy_size = chunk->num_columns;
    data->z_size = 2*chunk->num_z_cells;
  }
  else if (centering == PRISMESH_XYEDGE)
  {
    data->xy_size = 2*chunk->num_xy_faces;
    data->z_size = 2*chunk->num_z_cells;
  }
  else if (centering == PRISMESH_ZEDGE)
  {
    data->xy_size = 2*chunk->num_xy_faces;
    data->z_size = 2*chunk->num_z_cells;
  }
  else // (centering == PRISMESH_NODE)
  {
    data->xy_size = chunk->num_xy_faces;
    data->z_size = 2*chunk->num_z_cells;
  }
  size_t data_size = data->xy_size * data->z_size * num_components;
  data->data = polymec_calloc(sizeof(real_t) * data_size);
  return data;
}

static void prismesh_chunk_data_free(prismesh_chunk_data_t* data)
{
  polymec_free(data->data);
  polymec_free(data);
}

void prismesh_chunk_data_copy(prismesh_chunk_data_t* data, 
                              prismesh_chunk_data_t* dest)
{
  ASSERT(dest->chunk == data->chunk);
  ASSERT(dest->centering == data->centering);
  ASSERT(dest->num_components == data->num_components);
  ASSERT(dest->xy_size == data->xy_size);
  ASSERT(dest->z_size == data->z_size);
  memcpy(dest->data, data->data, prismesh_chunk_data_size(data));
}

DEFINE_UNORDERED_MAP(chunk_data_map, int, prismesh_chunk_data_t*, int_hash, int_equals)

struct prismesh_field_t 
{
  prismesh_t* mesh;
  prismesh_centering_t centering;
  size_t num_components;

  size_t num_xy_chunks, num_z_chunks;
  chunk_data_map_t* chunks;
};

static inline int chunk_index(prismesh_field_t* field, int xy_index, int z_index)
{
  return (int)(field->num_z_chunks * xy_index + z_index);
}

prismesh_field_t* prismesh_field_new(prismesh_t* mesh,
                                     prismesh_centering_t centering,
                                     size_t num_components)
{
  ASSERT(num_components > 0);
  prismesh_field_t* field = polymec_malloc(sizeof(prismesh_field_t));
  field->mesh = mesh;
  field->centering = centering;
  field->num_components = num_components;
  size_t num_xy_chunks, num_z_chunks, nz_per_chunk;
  prismesh_get_chunk_info(mesh, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  field->num_xy_chunks = num_xy_chunks;
  field->num_z_chunks = num_z_chunks;

  // Create data for each of the chunks in the mesh.
  int pos = 0, xy_index, z_index;
  prismesh_chunk_t* chunk;
  while (prismesh_next_chunk(mesh, &pos, &xy_index, &z_index, &chunk))
  {
    int ch_index = chunk_index(field, xy_index, z_index);
    prismesh_chunk_data_t* data = prismesh_chunk_data_new(chunk, centering, num_components);
    chunk_data_map_insert_with_v_dtor(field->chunks, ch_index, data, 
                                      prismesh_chunk_data_free);
  }

  return field;
}

void prismesh_field_free(prismesh_field_t* field)
{
  chunk_data_map_free(field->chunks);
  polymec_free(field);
}

prismesh_centering_t prismesh_field_centering(prismesh_field_t* field)
{
  return field->centering;
}

size_t prismesh_field_num_components(prismesh_field_t* field)
{
  return field->num_components;
}

size_t prismesh_field_num_chunks(prismesh_field_t* field)
{
  return field->chunks->size;
}

prismesh_t* prismesh_field_mesh(prismesh_field_t* field)
{
  return field->mesh;
}

void prismesh_field_copy(prismesh_field_t* field,
                         prismesh_field_t* dest)
{
  ASSERT(dest->mesh == field->mesh);
  ASSERT(dest->centering == field->centering);
  ASSERT(dest->num_components == field->num_components);
  ASSERT(dest->num_xy_chunks == field->num_xy_chunks);
  ASSERT(dest->num_z_chunks == field->num_z_chunks);
  int pos = 0, ch_index;
  prismesh_chunk_data_t* src_data;
  while (chunk_data_map_next(field->chunks, &pos, &ch_index, &src_data))
  {
    prismesh_chunk_data_t** dest_data_p = chunk_data_map_get(dest->chunks, ch_index);
    ASSERT(dest_data_p != NULL);
    prismesh_chunk_data_t* dest_data = *dest_data_p;
    ASSERT(dest_data->chunk == src_data->chunk);
    prismesh_chunk_data_copy(src_data, dest_data);
  }
}

prismesh_chunk_data_t* prismesh_field_chunk_data(prismesh_field_t* field, 
                                                 int xy_index,
                                                 int z_index)
{
  int index = chunk_index(field, xy_index, z_index);
  prismesh_chunk_data_t** data_p = chunk_data_map_get(field->chunks, index);
  if (data_p != NULL)
    return *data_p;
  else
    return NULL;
}

bool prismesh_field_next_chunk(prismesh_field_t* field, int* pos, 
                               int* xy_index, int* z_index,
                               prismesh_chunk_data_t** chunk_data)
{
  prismesh_chunk_t* chunk;
  bool result = prismesh_next_chunk(field->mesh, pos, xy_index, z_index, &chunk);
  if (result)
  {
    *chunk_data = prismesh_field_chunk_data(field, *xy_index, *z_index);
    ASSERT((*chunk_data)->chunk == chunk);
  }
  return result;
}

real_enumerable_generator_t* prismesh_field_enumerate(prismesh_field_t* field)
{
  size_t num_values = 0;
  int pos = 0, xy, z;
  prismesh_chunk_data_t* chunk_data;
  while (prismesh_field_next_chunk(field, &pos, &xy, &z, &chunk_data))
    num_values += prismesh_chunk_data_size(chunk_data) / sizeof(real_t);

  size_t offset = 0;
  real_t* array = polymec_malloc(sizeof(real_t) * num_values);
  while (prismesh_field_next_chunk(field, &pos, &xy, &z, &chunk_data))
  {
    size_t data_size = prismesh_chunk_data_size(chunk_data);
    memcpy(&array[offset], chunk_data->data, data_size);
    offset += data_size / sizeof(real_t);
  }
  return real_enumerable_generator_from_array(array, num_values, true);
}

real_enumerable_generator_t* prismesh_chunk_data_enumerate(prismesh_chunk_data_t* chunk_data)
{
  size_t num_values = prismesh_chunk_data_size(chunk_data) / sizeof(real_t);
  return real_enumerable_generator_from_array((real_t*)chunk_data->data, num_values, NULL);
}

