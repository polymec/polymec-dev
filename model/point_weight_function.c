// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/point_weight_function.h"

struct point_weight_function_t 
{
  char* name;
  void* context;
  point_weight_function_vtable vtable;
};

point_weight_function_t* point_weight_function_new(const char* name,
                                                   void* context,
                                                   point_weight_function_vtable vtable)
{
  ASSERT(vtable.value != NULL);
  ASSERT(vtable.gradient != NULL);
  point_weight_function_t* W = polymec_malloc(sizeof(point_weight_function_t));
  W->name = string_dup(name);
  W->context = context;
  W->vtable = vtable;
  return W;
}

void point_weight_function_free(point_weight_function_t* W)
{
  if ((W->vtable.dtor != NULL) && (W->context != NULL))
    W->vtable.dtor(W->context);
  string_free(W->name);
  polymec_free(W);
}

const char* point_weight_function_name(point_weight_function_t* W)
{
  return (const char*)W->name;
}

real_t point_weight_function_value(point_weight_function_t* W, vector_t* y)
{
  return W->vtable.value(W->context, y);
}

// Returns the gradient of the point weight function for the given displacement 
// vector y = x - x0.
vector_t point_weight_function_gradient(point_weight_function_t* W, vector_t* y)
{
  return W->vtable.gradient(W->context, y);
}

