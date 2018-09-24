// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/point_cloud_field.h"

// Observer function.
static void pcf_set_num_ghosts(void* context, size_t num_ghosts)
{
  point_cloud_field_t* field = context;
  field->num_ghost_values = num_ghosts;
  if (field->num_local_values + field->num_ghost_values > field->capacity)
  {
    field->capacity = field->num_local_values + field->num_ghost_values;
    field->data = polymec_realloc(field->data, sizeof(real_t) * field->num_components * field->capacity);
  }
}

point_cloud_field_t* point_cloud_field_new(point_cloud_t* cloud,
                                           size_t num_components)
{
  ASSERT(num_components > 0);
  point_cloud_field_t* field = polymec_malloc(sizeof(point_cloud_field_t));
  field->cloud = cloud;
  field->num_components = num_components;
  field->num_local_values= cloud->num_points;
  field->num_ghost_values = cloud->num_ghosts;
  field->capacity = field->num_local_values + field->num_ghost_values;
  field->data = polymec_calloc(sizeof(real_t) * num_components * field->capacity);

  // Register this field as an observer on its cloud.
  point_cloud_observer_vtable obs_vtable = {.set_num_ghosts = pcf_set_num_ghosts};
  point_cloud_observer_t* obs = point_cloud_observer_new(field, obs_vtable);
  point_cloud_add_observer(cloud, obs);

  return field;
}

void point_cloud_field_free(point_cloud_field_t* field)
{
  polymec_free(field->data);
  polymec_free(field);
}

bool point_cloud_field_compare_all(point_cloud_field_t* field,
                                   point_cloud_field_t* other_field,
                                   int component,
                                   bool (*comparator)(real_t val, real_t other_val))
{
  // We need a fair comparison.
  ASSERT(field->num_local_values == other_field->num_local_values);
  ASSERT(field->num_components > (size_t)component);
  ASSERT(other_field->num_components > (size_t)component);

  bool all = true;
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f, field);
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f1, other_field);
  for (size_t i = 0; i < field->num_local_values; ++i)
  {
    all = all && comparator(f[i][component], f1[i][component]);
    if (!all) break;
  }
  return all;
}

bool point_cloud_field_compare_any(point_cloud_field_t* field,
                                   point_cloud_field_t* other_field,
                                   int component,
                                   bool (*comparator)(real_t val, real_t other_val))
{
  // We need a fair comparison.
  ASSERT(field->num_local_values == other_field->num_local_values);
  ASSERT(field->num_components > (size_t)component);
  ASSERT(other_field->num_components > (size_t)component);

  bool any = false;
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f, field);
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f1, other_field);
  for (size_t i = 0; i < field->num_local_values; ++i)
  {
    any = comparator(f[i][component], f1[i][component]);
    if (any) break;
  }
  return any;
}

bool point_cloud_field_compare_none(point_cloud_field_t* field,
                                    point_cloud_field_t* other_field,
                                    int component,
                                    bool (*comparator)(real_t val, real_t other_val))
{
  // We need a fair comparison.
  ASSERT(field->num_local_values == other_field->num_local_values);
  ASSERT(field->num_components > (size_t)component);
  ASSERT(other_field->num_components > (size_t)component);

  bool none = true;
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f, field);
  DECLARE_POINT_CLOUD_FIELD_ARRAY(f1, other_field);
  for (size_t i = 0; i < field->num_local_values; ++i)
  {
    none = !comparator(f[i][component], f1[i][component]);
    if (!none) break;
  }
  return none;
}

