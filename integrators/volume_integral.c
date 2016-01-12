// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "integrators/volume_integral.h"

struct volume_integral_t 
{
  char* name;
  void* context;
  volume_integral_vtable vtable;

  // Quadrature points/weights.
  int num_points;
  point_t* points;
  real_t* weights;
};

static void volume_integral_free(void* ctx, void* dummy)
{
  volume_integral_t* integ = ctx;
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  if (integ->points != NULL)
    polymec_free(integ->points);
  if (integ->weights != NULL)
    polymec_free(integ->weights);
  string_free(integ->name);
}

volume_integral_t* volume_integral_new(const char* name,
                                       void* context,
                                       volume_integral_vtable vtable)
{
  ASSERT(vtable.num_quad_points != NULL);
  ASSERT(vtable.get_quadrature != NULL);

  volume_integral_t* integ = GC_MALLOC(sizeof(volume_integral_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->num_points = 0;
  integ->points = NULL;
  integ->weights = NULL;
  GC_register_finalizer(integ, volume_integral_free, integ, NULL, NULL);
  return integ;
}

void volume_integral_set_domain(volume_integral_t* integ, int i)
{
  integ->num_points = integ->vtable.num_quad_points(integ->context, i);
  if (integ->num_points > 0)
  {
    integ->points = polymec_realloc(integ->points, sizeof(point_t) * integ->num_points);
    integ->weights = polymec_realloc(integ->weights, sizeof(real_t) * integ->num_points);
    integ->vtable.get_quadrature(integ->context, i, integ->points, integ->weights);
  }
}

void volume_integral_compute(volume_integral_t* integ, sp_func_t* f, real_t* integral)
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

void volume_integral_compute_at_time(volume_integral_t* integ, real_t t, st_func_t* f, real_t* integral)
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

int volume_integral_num_points(volume_integral_t* integ)
{
  return integ->num_points;
}

void volume_integral_get_quadrature(volume_integral_t* integ,
                                    point_t* points,
                                    real_t* weights)
{
  if (integ->num_points > 0)
  {
    memcpy(points, integ->points, sizeof(point_t) * integ->num_points);
    memcpy(weights, integ->weights, sizeof(real_t) * integ->num_points);
  }
}


bool volume_integral_next_quad_point(volume_integral_t* integ, int* pos,
                                     point_t* point, real_t* weight)
{
  if (*pos >= integ->num_points)
    return false;
  *point = integ->points[*pos];
  *weight = integ->weights[*pos];
  ++(*pos);
  return true;
}

