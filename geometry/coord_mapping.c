// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"
#include "geometry/coord_mapping.h"

struct coord_mapping_t 
{
  char* name;
  void* context;
  coord_mapping_vtable vtable;
};

static void coord_mapping_free(void* ctx)
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
  ASSERT(vtable.jacobian != NULL);
  coord_mapping_t* m = polymec_gc_malloc(sizeof(coord_mapping_t), coord_mapping_free);
  m->name = string_dup(name);
  m->context = context;
  m->vtable = vtable;
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
  if (mapping->vtable.map_vector != NULL)
    mapping->vtable.map_vector(mapping->context, x, v, v1);
  else
  {
    // v1 = J * v.
    tensor2_t J;
    mapping->vtable.jacobian(mapping->context, x, &J);
    tensor2_dot_vector(&J, v, v1);
  }
}

void coord_mapping_compute_jacobian(coord_mapping_t* mapping, point_t* x, tensor2_t* J)
{
  mapping->vtable.jacobian(mapping->context, x, J);
}

coord_mapping_t* coord_mapping_inverse(coord_mapping_t* mapping)
{
  if (mapping->vtable.inverse != NULL)
    return mapping->vtable.inverse(mapping->context);
  else
    return NULL;
}

real_t coord_mapping_det_J(coord_mapping_t* mapping, point_t* x)
{
  if (mapping->vtable.det_J != NULL)
    return mapping->vtable.det_J(mapping->context, x);
  else
  {
    tensor2_t J;
    mapping->vtable.jacobian(mapping->context, x, &J);
    return tensor2_det(&J);
  }
}

void coord_mapping_compute_metric(coord_mapping_t* mapping, point_t* x, tensor2_t* G)
{
  if (mapping->vtable.metric != NULL)
    mapping->vtable.metric(mapping->context, x, G);
  else
  {
    // The metric G is just J^T * J.
    tensor2_t J;
    mapping->vtable.jacobian(mapping->context, x, &J);
    char trans = 'T', no_trans = 'N';
    int three = 3;
    real_t one = 1.0, zero = 0.0;
    rgemm(&trans, &no_trans, &three, &three, &three, &one, (real_t*)&J, 
          &three, (real_t*)&J, &three, &zero, (real_t*)G, &three);
  }
}

typedef struct
{
  coord_mapping_t* map1;
  coord_mapping_t* map2;
} comp_cm_t;

static void comp_map_point(void* context, point_t* x, point_t* y)
{
  comp_cm_t* comp = context;
  point_t x1;
  coord_mapping_map_point(comp->map1, x, &x1);
  coord_mapping_map_point(comp->map2, &x1, y);
}

static void comp_jacobian(void* context, point_t* x, tensor2_t* J)
{
  // J = J1 * J2.
  comp_cm_t* comp = context;
  tensor2_t J1, J2;
  coord_mapping_compute_jacobian(comp->map1, x, &J1);
  coord_mapping_compute_jacobian(comp->map1, x, &J2);
  char no_trans = 'N';
  int three = 3;
  real_t one = 1.0, zero = 0.0;
  rgemm(&no_trans, &no_trans, &three, &three, &three, &one, (real_t*)&J1, 
        &three, (real_t*)&J2, &three, &zero, (real_t*)J, &three);
}

coord_mapping_t* composite_coord_mapping_new(coord_mapping_t* map1, coord_mapping_t* map2)
{
  comp_cm_t* comp = polymec_malloc(sizeof(comp_cm_t));
  comp->map1 = map1;
  comp->map2 = map2;
  char name[2048];
  snprintf(name, 2047, "%s o %s", map1->name, map2->name);
  coord_mapping_vtable vtable = {.map_point = comp_map_point,
                                 .jacobian = comp_jacobian,
                                 .dtor = polymec_free};
  return coord_mapping_new(name, comp, vtable);
}

