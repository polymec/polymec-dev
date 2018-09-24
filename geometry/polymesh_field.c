// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/polymesh_field.h"

// Constructs a new polymesh field with the given number of components
// on the given mesh.
polymesh_field_t* polymesh_field_new(polymesh_t* mesh,
                                     polymesh_centering_t centering,
                                     size_t num_components)
{
  ASSERT(num_components > 0);
  polymesh_field_t* field = polymec_malloc(sizeof(polymesh_field_t));
  field->mesh = mesh;
  field->centering = centering;
  field->num_components = num_components;
  switch (centering)
  {
    case POLYMESH_CELL: field->num_local_values = mesh->num_cells; break;
    case POLYMESH_FACE: field->num_local_values = mesh->num_faces; break;
    case POLYMESH_EDGE: field->num_local_values = mesh->num_edges; break;
    case POLYMESH_NODE: field->num_local_values = mesh->num_nodes;
  }
  field->num_ghost_values = (centering == POLYMESH_CELL) ? mesh->num_ghost_cells : 0;
  field->capacity = field->num_local_values + field->num_ghost_values;
  field->data = polymec_calloc(sizeof(real_t) * num_components * field->capacity);

  return field;
}

// Destroys the given polymesh field.
void polymesh_field_free(polymesh_field_t* field)
{
  polymec_free(field->data);
  polymec_free(field);
}

bool polymesh_field_compare_all(polymesh_field_t* field,
                                polymesh_field_t* other_field,
                                int component,
                                bool (*comparator)(real_t val, real_t other_val))
{
  // We need a fair comparison.
  ASSERT(field->mesh == other_field->mesh);
  ASSERT(field->centering == other_field->centering);
  ASSERT(field->num_components > (size_t)component);
  ASSERT(other_field->num_components > (size_t)component);

  bool all = true;
  DECLARE_POLYMESH_FIELD_ARRAY(f, field);
  DECLARE_POLYMESH_FIELD_ARRAY(f1, other_field);
  for (size_t i = 0; i < field->num_local_values; ++i)
  {
    all = all && comparator(f[i][component], f1[i][component]);
    if (!all) break;
  }
  return all;
}

bool polymesh_field_compare_any(polymesh_field_t* field,
                                polymesh_field_t* other_field,
                                int component,
                                bool (*comparator)(real_t val, real_t other_val))
{
  // We need a fair comparison.
  ASSERT(field->mesh == other_field->mesh);
  ASSERT(field->centering == other_field->centering);
  ASSERT(field->num_components > (size_t)component);
  ASSERT(other_field->num_components > (size_t)component);

  bool any = false;
  DECLARE_POLYMESH_FIELD_ARRAY(f, field);
  DECLARE_POLYMESH_FIELD_ARRAY(f1, other_field);
  for (size_t i = 0; i < field->num_local_values; ++i)
  {
    any = comparator(f[i][component], f1[i][component]);
    if (any) break;
  }
  return any;
}

bool polymesh_field_compare_none(polymesh_field_t* field,
                                 polymesh_field_t* other_field,
                                 int component,
                                 bool (*comparator)(real_t val, real_t other_val))
{
  // We need a fair comparison.
  ASSERT(field->mesh == other_field->mesh);
  ASSERT(field->centering == other_field->centering);
  ASSERT(field->num_components > (size_t)component);
  ASSERT(other_field->num_components > (size_t)component);

  bool none = true;
  DECLARE_POLYMESH_FIELD_ARRAY(f, field);
  DECLARE_POLYMESH_FIELD_ARRAY(f1, other_field);
  for (size_t i = 0; i < field->num_local_values; ++i)
  {
    none = !comparator(f[i][component], f1[i][component]);
    if (!none) break;
  }
  return none;
}

