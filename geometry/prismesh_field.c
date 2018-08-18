// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/prismesh_field.h"

static prismesh_chunk_data_t* prismesh_chunk_data_new(prismesh_chunk_t* chunk,
                                                      prismesh_centering_t centering,
                                                      size_t num_components)
{
  prismesh_chunk_data_t* data = polymec_malloc(sizeof(prismesh_chunk_data_t));
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
  data->data = polymec_malloc(sizeof(real_t) * data_size);
  memset(data->data, 0, sizeof(real_t) * data_size);
  return data;
}

static void prismesh_chunk_data_free(prismesh_chunk_data_t* data)
{
  polymec_free(data->data);
  polymec_free(data);
}

DEFINE_ARRAY(chunk_data_array, prismesh_chunk_data_t*)

struct prismesh_field_t 
{
  prismesh_t* mesh;
  prismesh_centering_t centering;
  size_t num_components;

  int_array_t* chunks;
  chunk_data_array_t* chunk_data;
};

prismesh_field_t* prismesh_field_new(prismesh_t* mesh,
                                     prismesh_centering_t centering,
                                     size_t num_components)
{
  ASSERT(num_components > 0);
  prismesh_field_t* field = polymec_malloc(sizeof(prismesh_field_t));
  field->mesh = mesh;
  field->centering = centering;
  field->num_components = num_components;

  // Create data for each of the chunks in the mesh.
  int pos = 0;
  prismesh_chunk_t* chunk;
  while (prismesh_next_chunk(mesh, &pos, &chunk))
  {
    prismesh_chunk_data_t* data = prismesh_chunk_data_new(chunk, centering, num_components);
    chunk_data_array_append_with_dtor(field->chunk_data, data, prismesh_chunk_data_free);
  }

  return field;
}

void prismesh_field_free(prismesh_field_t* field)
{
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

bool prismesh_field_next_chunk(prismesh_field_t* field, int* pos, 
                               prismesh_chunk_t** chunk, 
                               prismesh_chunk_data_t** data)
{
  if (*pos >= (int)field->chunks->size)
    return false;
  else
  {
    *data = field->chunk_data->data[*pos];
    return prismesh_next_chunk(field->mesh, pos, chunk);
  }
}

