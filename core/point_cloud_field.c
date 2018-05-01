// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/point_cloud_field.h"

struct point_cloud_field_t 
{
  point_cloud_t* cloud;
  int num_comps, num_local, num_ghost, capacity;
  real_t* data;
};

// Observer function.
static void pcf_set_num_ghosts(void* context, int num_ghosts)
{
  point_cloud_field_t* field = context;
  field->num_ghost = num_ghosts;
  if (field->num_local + field->num_ghost > field->capacity)
  {
    field->capacity = field->num_local + field->num_ghost;
    field->data = polymec_realloc(field->data, sizeof(real_t) * field->num_comps * field->capacity);
  }
}

point_cloud_field_t* point_cloud_field_new(point_cloud_t* cloud,
                                           int num_components)
{
  ASSERT(num_components > 0);
  point_cloud_field_t* field = polymec_malloc(sizeof(point_cloud_field_t));
  field->cloud = cloud;
  field->num_comps = num_components;
  field->num_local = cloud->num_points;
  field->num_ghost = cloud->num_ghosts;
  field->capacity = field->num_local + field->num_ghost;
  int N = field->num_local + field->num_ghost;
  field->data = polymec_calloc(sizeof(real_t) * num_components * N);

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

point_cloud_t* point_cloud_field_cloud(point_cloud_field_t* field)
{
  return field->cloud;
}

int point_cloud_field_num_components(point_cloud_field_t* field)
{
  return field->num_comps;
}

int point_cloud_field_num_local_values(point_cloud_field_t* field)
{
  return field->num_local;
}

int point_cloud_field_num_ghost_values(point_cloud_field_t* field)
{
  return field->num_ghost;
}

real_t* point_floud_field_data(point_cloud_field_t* field)
{
  return field->data;
}

