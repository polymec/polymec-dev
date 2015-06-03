// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/shape_function.h"

struct kernel_function_t
{
  char* name;
  void* context;
  void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* value, vector_t* gradient);
  void (*dtor)(void* context);
};

kernel_function_t* kernel_function_new(const char* name,
                                       void* context,
                                       void (*compute)(void* context, point_t* points, real_t* extents, int num_points, point_t* x, real_t* value, vector_t* gradient),
                                       void (*dtor)(void* context))
{
  kernel_function_t* W = polymec_malloc(sizeof(kernel_function_t));
  W->name = string_dup(name);
  W->context = context;
  W->compute = compute;
  W->dtor = dtor;
  return W;
}

void kernel_function_free(kernel_function_t* W)
{
  string_free(W->name);
  if ((W->dtor != NULL) && (W->context != NULL))
    W->dtor(W->context);
  polymec_free(W);
}

void kernel_function_compute(kernel_function_t* W, 
                             point_t* points, 
                             real_t* extents, 
                             int num_points, 
                             point_t* x, 
                             real_t* values,
                             vector_t* gradients)
{
  W->compute(W->context, points, extents, num_points, x, values, gradients);
}

struct shape_function_t 
{
  char* name;
  void* context;
  shape_function_vtable vtable;
  int i, N;
};

shape_function_t* shape_function_new(const char* name, 
                                     void* context, 
                                     shape_function_vtable vtable)
{
  ASSERT(vtable.neighborhood_size != NULL);
  ASSERT(vtable.get_neighborhood_points != NULL);
  ASSERT(vtable.compute != NULL);

  shape_function_t* phi = polymec_malloc(sizeof(shape_function_t));
  phi->name = string_dup(name);
  phi->context = context;
  phi->vtable = vtable;
  phi->i = -1;
  phi->N = -1;
  return phi;
}

void shape_function_free(shape_function_t* phi)
{
  string_free(phi->name);
  if ((phi->vtable.dtor != NULL) && (phi->context != NULL))
    phi->vtable.dtor(phi->context);
  polymec_free(phi);
}

void shape_function_set_neighborhood(shape_function_t* phi, int point_index)
{
  phi->i = point_index;
#if 0
  point_t* xi = &phi->domain->points[phi->i];
  int Nj = stencil_size(phi->neighborhoods, phi->i);
  point_t xj[Nj];
  int pos = 0, j, k = 0;
  while (stencil_next(phi->neighborhoods, phi->i, &pos, &j, NULL))
    xj[k++] = phi->domain->points[
#endif
  phi->N = phi->vtable.neighborhood_size(phi->context, point_index);
  if (phi->vtable.set_neighborhood != NULL)
    phi->vtable.set_neighborhood(phi->context, phi->i);
}

int shape_function_num_points(shape_function_t* phi)
{
  return phi->N;
}

void shape_function_get_points(shape_function_t* phi, point_t* points)
{
  phi->vtable.get_neighborhood_points(phi->context, phi->i, points);
}

void shape_function_compute(shape_function_t* phi, 
                            point_t* x,
                            real_t* values,
                            vector_t* gradients)
{
  phi->vtable.compute(phi->context, phi->i, x, values, gradients);
}

