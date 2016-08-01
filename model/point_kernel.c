// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "model/point_kernel.h"

struct point_kernel_t
{
  char* name;
  void* context;
  void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* value, vector_t* gradient);
  void (*dtor)(void* context);
};

static void point_kernel_free(void* ctx, void* dummy)
{
  point_kernel_t* kernel = ctx;
  string_free(kernel->name);
  if ((kernel->dtor != NULL) && (kernel->context != NULL))
    kernel->dtor(kernel->context);
}

point_kernel_t* point_kernel_new(const char* name,
                                 void* context,
                                 void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* values, vector_t* gradients),
                                 void (*dtor)(void* context))
{
  point_kernel_t* kernel = GC_MALLOC(sizeof(point_kernel_t));
  kernel->name = string_dup(name);
  kernel->context = context;
  kernel->compute = compute;
  kernel->dtor = dtor;
  GC_register_finalizer(kernel, point_kernel_free, kernel, NULL, NULL);
  return kernel;
}

void point_kernel_compute(point_kernel_t* kernel, 
                          point_t* points, 
                          real_t* extents, 
                          int num_points, 
                          point_t* x, 
                          real_t* values,
                          vector_t* gradients)
{
  kernel->compute(kernel->context, points, extents, num_points, x, 
                  values, gradients);
}

static void spline4_compute(void* context, 
                            point_t* points, 
                            real_t* extents, 
                            int num_points, 
                            point_t* x, 
                            real_t* values, 
                            vector_t* gradients)
{
  real_t eta_max = *((real_t*)context);
  for (int i = 0; i < num_points; ++i)
  {
    point_t* xi = &points[i];
    real_t h = extents[i];
    real_t Di = point_distance(xi, x);
    real_t eta = Di/h;
    real_t X = eta/eta_max;
    if (X < 1.0)
    {
      values[i] = 1.0 - 6.0*X*X + 8.0*X*X*X - 3.0*X*X*X*X;
      real_t dWdDi = (-12.0*X + 24.0*X*X - 12.0*X*X*X)/(eta_max*h);
      gradients[i].x = dWdDi * (xi->x-x->x)/Di;
      gradients[i].y = dWdDi * (xi->y-x->y)/Di;
      gradients[i].z = dWdDi * (xi->z-x->z)/Di;
    }
    else
    {
      values[i] = 0.0;
      gradients[i].x = 0.0;
      gradients[i].y = 0.0;
      gradients[i].z = 0.0;
    }
  }
}

point_kernel_t* spline4_point_kernel_new(real_t eta_max)
{
  char name[1024];
  snprintf(name, 1023, "Spline4 kernel");
  real_t* context = polymec_malloc(sizeof(real_t));
  *context = eta_max;
  return point_kernel_new((char*)name, context, spline4_compute, polymec_free);
}

