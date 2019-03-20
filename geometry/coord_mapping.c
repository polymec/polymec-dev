// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"
#include "geometry/coord_mapping.h"
#include "geometry/field_metadata.h"

struct coord_mapping_t
{
  char* name;
  void* context;
  coord_mapping_vtable vtable;
  coord_mapping_t* inverse;
};

static void coord_mapping_free(void* ctx)
{
  coord_mapping_t* mapping = (coord_mapping_t*)ctx;
  if (mapping->inverse != NULL)
    release_ref(mapping->inverse);
  if (mapping->vtable.dtor != NULL)
    free(mapping->context);
  free(mapping->name);
}

coord_mapping_t* coord_mapping_new(const char* name,
                                   void* context,
                                   coord_mapping_vtable vtable)
{
  ASSERT(context != NULL);
  ASSERT(vtable.map_point != NULL);
  ASSERT(vtable.jacobian != NULL);
  coord_mapping_t* m = polymec_refcounted_malloc(sizeof(coord_mapping_t),
                                                 coord_mapping_free);
  m->name = string_dup(name);
  m->context = context;
  m->vtable = vtable;
  m->inverse = NULL;
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

void coord_mapping_map_point(coord_mapping_t* mapping, point_t* x,
                             point_t* y)
{
  mapping->vtable.map_point(mapping->context, x, y);
}

void coord_mapping_map_vector(coord_mapping_t* mapping, point_t* x,
                              vector_t* v, vector_t* v1)
{
  if (mapping->vtable.map_vector != NULL)
    mapping->vtable.map_vector(mapping->context, x, v, v1);
  else
  {
    // v1 = J * v.
    real_t J[3][3];
    mapping->vtable.jacobian(mapping->context, x, J);
    v1->x = J[0][0]*v->x + J[0][1]*v->y + J[0][2]*v->z;
    v1->y = J[1][0]*v->x + J[1][1]*v->y + J[1][2]*v->z;
    v1->z = J[2][0]*v->x + J[2][1]*v->y + J[2][2]*v->z;
  }
}

void coord_mapping_map_tensor2(coord_mapping_t* mapping, point_t* x,
                               tensor2_t* t, tensor2_t* t1)
{
  if (mapping->vtable.map_tensor2 != NULL)
    mapping->vtable.map_tensor2(mapping->context, x, t, t1);
  else
  {
    real_t J[3][3];
    mapping->vtable.jacobian(mapping->context, x, J);
    t1->xx = J[0][0]*(J[0][0]*t->xx + J[0][1]*t->xy + J[0][2]*t->xz) +
             J[0][1]*(J[0][0]*t->yx + J[0][1]*t->yy + J[0][2]*t->yz) +
             J[0][2]*(J[0][0]*t->yx + J[0][1]*t->yy + J[0][2]*t->yz);
    t1->xy = J[0][0]*(J[1][0]*t->xx + J[1][1]*t->xy + J[1][2]*t->xz) +
             J[0][1]*(J[1][0]*t->yx + J[1][1]*t->yy + J[1][2]*t->yz) +
             J[0][2]*(J[1][0]*t->yx + J[1][1]*t->yy + J[1][2]*t->yz);
    t1->xz = J[0][0]*(J[2][0]*t->xx + J[2][1]*t->xy + J[2][2]*t->xz) +
             J[0][1]*(J[2][0]*t->yx + J[2][1]*t->yy + J[2][2]*t->yz) +
             J[0][2]*(J[2][0]*t->yx + J[2][1]*t->yy + J[2][2]*t->yz);

    t1->yx = J[1][0]*(J[0][0]*t->xx + J[0][1]*t->xy + J[0][2]*t->xz) +
             J[1][1]*(J[0][0]*t->yx + J[0][1]*t->yy + J[0][2]*t->yz) +
             J[1][2]*(J[0][0]*t->yx + J[0][1]*t->yy + J[0][2]*t->yz);
    t1->yy = J[1][0]*(J[1][0]*t->xx + J[1][1]*t->xy + J[1][2]*t->xz) +
             J[1][1]*(J[1][0]*t->yx + J[1][1]*t->yy + J[1][2]*t->yz) +
             J[1][2]*(J[1][0]*t->yx + J[1][1]*t->yy + J[1][2]*t->yz);
    t1->yz = J[1][0]*(J[2][0]*t->xx + J[2][1]*t->xy + J[2][2]*t->xz) +
             J[1][1]*(J[2][0]*t->yx + J[2][1]*t->yy + J[2][2]*t->yz) +
             J[1][2]*(J[2][0]*t->yx + J[2][1]*t->yy + J[2][2]*t->yz);

    t1->zx = J[2][0]*(J[0][0]*t->xx + J[0][1]*t->xy + J[0][2]*t->xz) +
             J[2][1]*(J[0][0]*t->yx + J[0][1]*t->yy + J[0][2]*t->yz) +
             J[2][2]*(J[0][0]*t->yx + J[0][1]*t->yy + J[0][2]*t->yz);
    t1->zy = J[2][0]*(J[1][0]*t->xx + J[1][1]*t->xy + J[1][2]*t->xz) +
             J[2][1]*(J[1][0]*t->yx + J[1][1]*t->yy + J[1][2]*t->yz) +
             J[2][2]*(J[1][0]*t->yx + J[1][1]*t->yy + J[1][2]*t->yz);
    t1->zz = J[2][0]*(J[2][0]*t->xx + J[2][1]*t->xy + J[2][2]*t->xz) +
             J[2][1]*(J[2][0]*t->yx + J[2][1]*t->yy + J[2][2]*t->yz) +
             J[2][2]*(J[2][0]*t->yx + J[2][1]*t->yy + J[2][2]*t->yz);
  }
}

void coord_mapping_map_symtensor2(coord_mapping_t* mapping, point_t* x,
                                  symtensor2_t* t, symtensor2_t* t1)
{
  if (mapping->vtable.map_symtensor2 != NULL)
    mapping->vtable.map_symtensor2(mapping->context, x, t, t1);
  else
  {
    real_t J[3][3];
    mapping->vtable.jacobian(mapping->context, x, J);
    t1->xx = J[0][0]*(J[0][0]*t->xx + J[0][1]*t->xy + J[0][2]*t->xz) +
             J[0][1]*(J[0][0]*t->xy + J[0][1]*t->yy + J[0][2]*t->yz) +
             J[0][2]*(J[0][0]*t->xy + J[0][1]*t->yy + J[0][2]*t->yz);
    t1->xy = J[0][0]*(J[1][0]*t->xx + J[1][1]*t->xy + J[1][2]*t->xz) +
             J[0][1]*(J[1][0]*t->xy + J[1][1]*t->yy + J[1][2]*t->yz) +
             J[0][2]*(J[1][0]*t->xy + J[1][1]*t->yy + J[1][2]*t->yz);
    t1->xz = J[0][0]*(J[2][0]*t->xx + J[2][1]*t->xy + J[2][2]*t->xz) +
             J[0][1]*(J[2][0]*t->xy + J[2][1]*t->yy + J[2][2]*t->yz) +
             J[0][2]*(J[2][0]*t->xy + J[2][1]*t->yy + J[2][2]*t->yz);
    t1->yy = J[1][0]*(J[1][0]*t->xx + J[1][1]*t->xy + J[1][2]*t->xz) +
             J[1][1]*(J[1][0]*t->xy + J[1][1]*t->yy + J[1][2]*t->yz) +
             J[1][2]*(J[1][0]*t->xy + J[1][1]*t->yy + J[1][2]*t->yz);
    t1->zz = J[2][0]*(J[2][0]*t->xx + J[2][1]*t->xy + J[2][2]*t->xz) +
             J[2][1]*(J[2][0]*t->xy + J[2][1]*t->yy + J[2][2]*t->yz) +
             J[2][2]*(J[2][0]*t->xy + J[2][1]*t->yy + J[2][2]*t->yz);
  }
}

void coord_mapping_compute_jacobian(coord_mapping_t* mapping,
                                    point_t* x,
                                    real_t J[3][3])
{
  mapping->vtable.jacobian(mapping->context, x, J);
}

coord_mapping_t* coord_mapping_inverse(coord_mapping_t* mapping)
{
  return mapping->inverse;
}

void coord_mapping_set_inverse(coord_mapping_t* mapping,
                               coord_mapping_t* inverse)
{
  ASSERT(inverse != NULL);
  mapping->inverse = inverse;
}

real_t coord_mapping_det_J(coord_mapping_t* mapping, point_t* x)
{
  if (mapping->vtable.det_J != NULL)
    return mapping->vtable.det_J(mapping->context, x);
  else
  {
    real_t J[3][3];
    mapping->vtable.jacobian(mapping->context, x, J);
    return J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) -
           J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) +
           J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[0][2]);
  }
}

void coord_mapping_compute_metric(coord_mapping_t* mapping, point_t* x,
                                  tensor2_t* G)
{
  if (mapping->vtable.metric != NULL)
    mapping->vtable.metric(mapping->context, x, G);
  else
  {
    // The metric G is just J^T * J.
    real_t J[3][3];
    mapping->vtable.jacobian(mapping->context, x, J);
    char trans = 'T', no_trans = 'N';
    int three = 3;
    real_t one = 1.0, zero = 0.0;
    rgemm(&trans, &no_trans, &three, &three, &three, &one, (real_t*)J,
          &three, (real_t*)J, &three, &zero, (real_t*)G, &three);
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

static void comp_jacobian(void* context, point_t* x, real_t J[3][3])
{
  // J = J1 * J2.
  comp_cm_t* comp = context;
  real_t J1[3][3], J2[3][3];
  coord_mapping_compute_jacobian(comp->map1, x, J1);
  coord_mapping_compute_jacobian(comp->map1, x, J2);
  char no_trans = 'N';
  int three = 3;
  real_t one = 1.0, zero = 0.0;
  rgemm(&no_trans, &no_trans, &three, &three, &three, &one, (real_t*)J1,
        &three, (real_t*)J2, &three, &zero, (real_t*)J, &three);
}

void coord_mapping_map_field_data(coord_mapping_t* mapping,
                                  field_metadata_t* metadata,
                                  point_t* x,
                                  real_t* field_data,
                                  real_t* mapped_field_data)
{
  int pos, c;

  // Map scalars.
  if (mapped_field_data != field_data)
  {
    pos = 0;
    while (field_metadata_next_scalar(metadata, &pos, &c))
      mapped_field_data[c] = field_data[c];
  }

  // Map vectors.
  pos = 0;
  while (field_metadata_next_vector(metadata, &pos, &c))
  {
    vector_t v = {field_data[c], field_data[c+1], field_data[c+2]}, v1;
    coord_mapping_map_vector(mapping, x, &v, &v1);
    mapped_field_data[c]   = v1.x;
    mapped_field_data[c+1] = v1.y;
    mapped_field_data[c+2] = v1.z;
  }

  // Map tensors.
  pos = 0;
  while (field_metadata_next_tensor2(metadata, &pos, &c))
  {
    tensor2_t t = {field_data[c],   field_data[c+1], field_data[c+2],
                   field_data[c+3], field_data[c+4], field_data[c+5],
                   field_data[c+6], field_data[c+7], field_data[c+8]}, t1;
    coord_mapping_map_tensor2(mapping, x, &t, &t1);
    mapped_field_data[c]   = t1.xx;
    mapped_field_data[c+1] = t1.xy;
    mapped_field_data[c+2] = t1.xz;
    mapped_field_data[c+3] = t1.yx;
    mapped_field_data[c+4] = t1.yy;
    mapped_field_data[c+5] = t1.yz;
    mapped_field_data[c+6] = t1.zx;
    mapped_field_data[c+7] = t1.zy;
    mapped_field_data[c+8] = t1.zz;
  }

  // Map symmetric tensors.
  pos = 0;
  while (field_metadata_next_symtensor2(metadata, &pos, &c))
  {
    symtensor2_t t = {field_data[c],   field_data[c+1], field_data[c+2],
                                       field_data[c+3], field_data[c+4],
                                                        field_data[c+5]}, t1;
    coord_mapping_map_symtensor2(mapping, x, &t, &t1);
    mapped_field_data[c]   = t1.xx;
    mapped_field_data[c+1] = t1.xy;
    mapped_field_data[c+2] = t1.xz;
    mapped_field_data[c+3] = t1.yy;
    mapped_field_data[c+4] = t1.yz;
    mapped_field_data[c+5] = t1.zz;
  }
}

coord_mapping_t* composite_coord_mapping_new(coord_mapping_t* map1,
                                             coord_mapping_t* map2)
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

