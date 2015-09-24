// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/point_spacing_estimator.h"

struct point_spacing_estimator_t 
{
  char* name;
  void* context;
  point_spacing_estimator_vtable vtable;
};

point_spacing_estimator_t* point_spacing_estimator_new(const char* name,
                                                       void* context,
                                                       point_spacing_estimator_vtable vtable)
{
  ASSERT(vtable.dx != NULL);
  point_spacing_estimator_t* estimator = polymec_malloc(sizeof(point_spacing_estimator_t));
  estimator->name = string_dup(name);
  estimator->context = context;
  estimator->vtable = vtable;
  return estimator;
}

void point_spacing_estimator_free(point_spacing_estimator_t* estimator)
{
  string_free(estimator->name);
  if ((estimator->context != NULL) && (estimator->vtable.dtor != NULL))
    estimator->vtable.dtor(estimator->context);
  polymec_free(estimator);
}

real_t point_spacing_estimator_dx(point_spacing_estimator_t* estimator, int i)
{
  return estimator->vtable.dx(estimator->context, i);
}

