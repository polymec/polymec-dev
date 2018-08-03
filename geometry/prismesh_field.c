// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/prismesh_field.h"

static prismesh_layer_data_t* prismesh_layer_data_new(prismesh_layer_t* layer,
                                                      prismesh_centering_t centering,
                                                      size_t num_components)
{
  prismesh_layer_data_t* data = polymec_malloc(sizeof(prismesh_layer_data_t));
  data->data = NULL;
  // FIXME
  return data;
}

static void prismesh_layer_data_free(prismesh_layer_data_t* data)
{
  polymec_free(data->data);
  polymec_free(data);
}

DEFINE_ARRAY(layer_data_array, prismesh_layer_data_t*)

struct prismesh_field_t 
{
  prismesh_t* mesh;
  prismesh_centering_t centering;
  size_t num_components;

  int_array_t* layers;
  layer_data_array_t* layer_data;
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

  // Create data for each of the layers in the mesh.
  int pos = 0;
  prismesh_layer_t* layer;
  while (prismesh_next_layer(mesh, &pos, &layer))
  {
    prismesh_layer_data_t* data = prismesh_layer_data_new(layer, centering, num_components);
    layer_data_array_append_with_dtor(field->layer_data, data, prismesh_layer_data_free);
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

size_t prismesh_field_num_layers(prismesh_field_t* field)
{
  return field->layers->size;
}

prismesh_t* prismesh_field_mesh(prismesh_field_t* field)
{
  return field->mesh;
}

bool prismesh_field_next_layer(prismesh_field_t* field, int* pos, 
                               prismesh_layer_t** layer, 
                               prismesh_layer_data_t** data)
{
  if (*pos >= (int)field->layers->size)
    return false;
  else
  {
    *data = field->layer_data->data[*pos];
    return prismesh_next_layer(field->mesh, pos, layer);
  }
}

