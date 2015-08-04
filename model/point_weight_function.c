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

vector_t point_weight_function_gradient(point_weight_function_t* W, vector_t* y)
{
  return W->vtable.gradient(W->context, y);
}

static real_t gaussian_value(void* context, vector_t* y)
{
  real_t epsilon = *((real_t*)context);
  real_t eps2 = epsilon*epsilon;
  real_t y2 = vector_dot(y, y);
  real_t exp_e2 = exp(-eps2);
  return (exp(-eps2*y2) - exp_e2) / (1.0 - exp_e2);
}

static vector_t gaussian_gradient(void* context, vector_t* y)
{
  real_t epsilon = *((real_t*)context);
  real_t eps2 = epsilon*epsilon;
  real_t y2 = vector_dot(y, y);
  real_t exp_e2 = exp(-eps2);
  real_t dWdu = -eps2 * exp(-eps2*y2) / (1.0 - exp_e2);
  vector_t grad = {.x = dWdu * y->x/y2, .y = dWdu * y->y/y2, .z = dWdu * y->z/y2};
  return grad;
}

point_weight_function_t* gaussian_point_weight_function_new(real_t epsilon)
{
  real_t* shape = polymec_malloc(sizeof(real_t));
  *shape = epsilon;
  char name[256];
  snprintf(name, 255, "Gaussian weight function (epsilon = %g)", epsilon);
  point_weight_function_vtable vtable = {.value = gaussian_value,
                                         .gradient = gaussian_gradient,
                                         .dtor = polymec_free};
  return point_weight_function_new(name, shape, vtable);
}
