// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/enumerable.h"
#include "core/timer.h"
#include "geometry/polymesh_field.h"

// Constructs a new polymesh field with the given number of components
// on the given mesh.
polymesh_field_t* polymesh_field_new(polymesh_t* mesh,
                                     polymesh_centering_t centering,
                                     int num_components)
{
  ASSERT(num_components > 0);
  polymesh_field_t* field = polymec_malloc(sizeof(polymesh_field_t));
  field->mesh = mesh;
  field->centering = centering;
  field->num_components = num_components;
  field->ex = NULL;
  switch (centering)
  {
    case POLYMESH_CELL: 
      field->num_local_values = mesh->num_cells; 
      break;
    case POLYMESH_FACE: 
      field->num_local_values = mesh->num_faces; 
      break;
    case POLYMESH_EDGE: 
      field->num_local_values = mesh->num_edges; 
      break;
    case POLYMESH_NODE: 
      field->num_local_values = mesh->num_nodes;
  }

  field->ex_token = -1;
  field->num_ghost_values = (centering == POLYMESH_CELL) ? mesh->num_ghost_cells : 0;
  field->capacity = field->num_local_values + field->num_ghost_values;
  field->data = polymec_calloc(num_components * field->capacity, sizeof(real_t));

  field->md = field_metadata_new(num_components);

  return field;
}

void polymesh_field_free(polymesh_field_t* field)
{
  if (field->ex != NULL)
    release_ref(field->ex);
  polymec_free(field->data);
  release_ref(field->md);
  polymec_free(field);
}

field_metadata_t* polymesh_field_metadata(polymesh_field_t* field)
{
  return field->md;
}

void polymesh_field_exchange(polymesh_field_t* field)
{
  polymesh_field_start_exchange(field);
  polymesh_field_finish_exchange(field);
}

void polymesh_field_start_exchange(polymesh_field_t* field)
{
  ASSERT(field->centering != POLYMESH_EDGE); // no edge exchanges!
  ASSERT(!polymesh_field_is_exchanging(field));
  START_FUNCTION_TIMER();

  // Do we have an exchanger yet?
  if (field->ex == NULL)
  {
    field->ex = polymesh_exchanger(field->mesh, field->centering);
    retain_ref(field->ex);
  }

  // Start the xy exchange.
  int stride = (int)field->num_components;
  field->ex_token = exchanger_start_exchange(field->ex, field->data, stride, 0, MPI_REAL_T);
  STOP_FUNCTION_TIMER();
}

void polymesh_field_finish_exchange(polymesh_field_t* field)
{
  ASSERT(polymesh_field_is_exchanging(field));
  START_FUNCTION_TIMER();
  if (field->ex_token != -1)
    exchanger_finish_exchange(field->ex, field->ex_token);
  field->ex_token = -1;
  STOP_FUNCTION_TIMER();
}

bool polymesh_field_is_exchanging(polymesh_field_t* field)
{
  return (field->ex_token != -1);
}

void polymesh_field_set_exchanger(polymesh_field_t* field, exchanger_t* ex)
{
  if (field->ex != NULL)
    release_ref(field->ex);
  if (ex != NULL)
    retain_ref(ex);
  field->ex = ex;
}

real_enumerable_generator_t* polymesh_field_enumerate(polymesh_field_t* field)
{
  size_t num_values = field->num_components * field->capacity;
  return real_enumerable_generator_from_array(field->data, num_values, NULL);
}

