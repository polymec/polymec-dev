// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "geometry/coord_mapping.h"

struct coord_mapping_t 
{
  char* name;
  void* context;
  coord_mapping_vtable vtable;
};

static void coord_mapping_free(void* ctx, void* dummy)
{
  coord_mapping_t* mapping = (coord_mapping_t*)ctx;
  if (mapping->vtable.dtor)
    free(mapping->context);
  free(mapping->name);
}

coord_mapping_t* coord_mapping_new(const char* name, void* context, coord_mapping_vtable vtable)
{
  ASSERT(context != NULL);
  ASSERT(vtable.map_point != NULL);
  ASSERT(vtable.map_vector != NULL);
  ASSERT(vtable.jacobian != NULL);
  coord_mapping_t* m = GC_MALLOC(sizeof(coord_mapping_t));
  m->name = string_dup(name);
  m->context = context;
  m->vtable = vtable;
  GC_register_finalizer(m, coord_mapping_free, m, NULL, NULL);
  return m;
}

const char* coord_mapping_name(coord_mapping_t* mapping)
{
  return mapping->name;
}

void* coord_mapping_context(coord_mapping_t* mapping)
{
  return mapping->context;
}

void coord_mapping_map_point(coord_mapping_t* mapping, point_t* x, point_t* y)
{
  mapping->vtable.map_point(mapping->context, x, y);
}

void coord_mapping_map_vector(coord_mapping_t* mapping, point_t* x, vector_t* v, vector_t* v1)
{
  mapping->vtable.map_vector(mapping->context, x, v, v1);
}

void coord_mapping_compute_jacobian(coord_mapping_t* mapping, point_t* x, real_t* J)
{
  mapping->vtable.jacobian(mapping->context, x, J);
}

