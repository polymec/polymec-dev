// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "model/point_spacing_estimator.h"

struct point_spacing_estimator_t 
{
  char* name;
  void* context;
  point_spacing_estimator_vtable vtable;
};

static void point_spacing_estimator_free(void* ctx, void* dummy)
{
  point_spacing_estimator_t* estimator = ctx;
  string_free(estimator->name);
  if ((estimator->context != NULL) && (estimator->vtable.dtor != NULL))
    estimator->vtable.dtor(estimator->context);
}

point_spacing_estimator_t* point_spacing_estimator_new(const char* name,
                                                       void* context,
                                                       point_spacing_estimator_vtable vtable)
{
  ASSERT(vtable.dx != NULL);
  point_spacing_estimator_t* estimator = GC_MALLOC(sizeof(point_spacing_estimator_t));
  estimator->name = string_dup(name);
  estimator->context = context;
  estimator->vtable = vtable;
  GC_register_finalizer(estimator, point_spacing_estimator_free, estimator, NULL, NULL);
  return estimator;
}

real_t point_spacing_estimator_dx(point_spacing_estimator_t* estimator, int i)
{
  return estimator->vtable.dx(estimator->context, i);
}

typedef struct
{
  point_cloud_t* points;
  stencil_t* stencil;
} stencil_dx_t;

static real_t stencil_dx(void* context, int i)
{
  stencil_dx_t* est = context;
  real_t dx = 0.0;
  int pos = 0, j, N = 0;
  point_t* xi = &(est->points->points[i]);
  while (stencil_next(est->stencil, i, &pos, &j, NULL))
  {
    point_t* xj = &(est->points->points[j]);
    dx += point_distance(xi, xj);
    ++N;
  }
  ASSERT(N > 0);
  dx /= N;
  return dx;
}

point_spacing_estimator_t* stencil_point_spacing_estimator_new(point_cloud_t* points,
                                                               stencil_t* stencil)
{
  ASSERT(points->num_points == stencil_num_indices(stencil));
  stencil_dx_t* est = polymec_malloc(sizeof(stencil_dx_t));
  est->points = points;
  est->stencil = stencil;
  point_spacing_estimator_vtable vtable = {.dx = stencil_dx, .dtor = polymec_free};
  return point_spacing_estimator_new("Stencil-based point spacing estimator",
                                     est, vtable);
}

