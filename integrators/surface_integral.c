// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "integrators/surface_integral.h"

struct surface_integral_t 
{
  char* name;
  void* context;
  surface_integral_vtable vtable;

  // Quadrature points/weights.
  int num_points;
  point_t* points;
  real_t* weights;
  vector_t* normals;
};

static void surface_integral_free(void* ctx, void* dummy)
{
  surface_integral_t* integ = ctx;
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  if (integ->points != NULL)
    polymec_free(integ->points);
  if (integ->weights != NULL)
    polymec_free(integ->weights);
  if (integ->normals != NULL)
    polymec_free(integ->normals);
  string_free(integ->name);
}

surface_integral_t* surface_integral_new(const char* name,
                                         void* context,
                                         surface_integral_vtable vtable)
{
  ASSERT(vtable.num_quad_points != NULL);
  ASSERT(vtable.get_quadrature != NULL);

  surface_integral_t* integ = GC_MALLOC(sizeof(surface_integral_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->num_points = 0;
  integ->points = NULL;
  integ->weights = NULL;
  integ->normals = NULL;
  GC_register_finalizer(integ, surface_integral_free, integ, NULL, NULL);
  return integ;
}

void surface_integral_set_domain(surface_integral_t* integ, int i)
{
  integ->num_points = integ->vtable.num_quad_points(integ->context, i);
  if (integ->num_points > 0)
  {
    integ->points = polymec_realloc(integ->points, sizeof(point_t) * integ->num_points);
    integ->weights = polymec_realloc(integ->weights, sizeof(real_t) * integ->num_points);
    integ->normals = polymec_realloc(integ->normals, sizeof(vector_t) * integ->num_points);
    integ->vtable.get_quadrature(integ->context, i, integ->points, integ->weights, integ->normals);
  }
}

void surface_integral_compute(surface_integral_t* integ, sp_func_t* f, real_t* integral)
{
  int num_comp = sp_func_num_comp(f);
  memset(integral, 0, num_comp * sizeof(real_t));
  for (int q = 0; q < integ->num_points; ++q)
  {
    real_t val[num_comp];
    sp_func_eval(f, &integ->points[q], val);
    for (int c = 0; c < num_comp; ++c)
      integral[c] += integ->weights[q] * val[c];
  }
}

void surface_integral_compute_at_time(surface_integral_t* integ, real_t t, st_func_t* f, real_t* integral)
{
  int num_comp = st_func_num_comp(f);
  memset(integral, 0, num_comp * sizeof(real_t));
  for (int q = 0; q < integ->num_points; ++q)
  {
    real_t val[num_comp];
    st_func_eval(f, &integ->points[q], t, val);
    for (int c = 0; c < num_comp; ++c)
      integral[c] += integ->weights[q] * val[c];
  }
}

real_t surface_integral_dot(surface_integral_t* integ, sp_func_t* f)
{
  int num_comp = sp_func_num_comp(f);
  ASSERT(num_comp == 3); // Must be a 3D vector function!
  real_t integral = 0.0;
  for (int q = 0; q < integ->num_points; ++q)
  {
    real_t val[num_comp];
    sp_func_eval(f, &integ->points[q], val);
    vector_t* n = &integ->normals[q];
    integral += integ->weights[q] * (val[0] * n->x + val[1] * n->y + val[2] * n->z);
  }
  return integral;
}

real_t surface_integral_dot_at_time(surface_integral_t* integ, 
                                    real_t t, st_func_t* f)
{
  int num_comp = st_func_num_comp(f);
  ASSERT(num_comp == 3); // Must be a 3D vector function!
  real_t integral = 0.0;
  for (int q = 0; q < integ->num_points; ++q)
  {
    real_t val[num_comp];
    st_func_eval(f, &integ->points[q], t, val);
    vector_t* n = &integ->normals[q];
    integral += integ->weights[q] * (val[0] * n->x + val[1] * n->y + val[2] * n->z);
  }
  return integral;
}

int surface_integral_num_points(surface_integral_t* integ)
{
  return integ->num_points;
}

void surface_integral_get_quadrature(surface_integral_t* integ,
                                     point_t* points,
                                     real_t* weights,
                                     vector_t* normals)
{
  if (integ->num_points > 0)
  {
    memcpy(points, integ->points, sizeof(point_t) * integ->num_points);
    memcpy(weights, integ->weights, sizeof(real_t) * integ->num_points);
    memcpy(normals, integ->normals, sizeof(vector_t) * integ->num_points);
  }
}


bool surface_integral_next_quad_point(surface_integral_t* integ, int* pos,
                                      point_t* point, real_t* weight, vector_t* normal)
{
  if (*pos >= integ->num_points)
    return false;
  *point = integ->points[*pos];
  *weight = integ->weights[*pos];
  *normal = integ->normals[*pos];
  ++(*pos);
  return true;
}

