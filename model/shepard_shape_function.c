// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/shepard_shape_function.h"

typedef struct
{
  int i;
  shape_function_kernel_t* W;
  point_cloud_t* domain;
  stencil_t* neighborhoods;
  real_t* smoothing_lengths;

  int N;
  point_t* xj;
  real_t* hj;
} shepard_t;

static int shepard_neighborhood_size(void* context, int i)
{
  shepard_t* shepard = context;
  return 1 + stencil_size(shepard->neighborhoods, i);
}

static void shepard_get_neighborhood_points(void* context, int i, point_t* points)
{
  shepard_t* shepard = context;
  int pos = 0, j, k = 0;
  points[k++] = shepard->domain->points[i];
  while (stencil_next(shepard->neighborhoods, i, &pos, &j, NULL))
    points[k++] = shepard->domain->points[j];
}

static void shepard_set_neighborhood(void* context, int i)
{
  shepard_t* shepard = context;
  ASSERT(i < shepard->domain->num_points); 

  // Extract the points.
  shepard->N = stencil_size(shepard->neighborhoods, i);
  int pos = 0, j, k = 0;
  while (stencil_next(shepard->neighborhoods, i, &pos, &j, NULL))
  {
    shepard->xj[k] = shepard->domain->points[j];
    shepard->hj[k] = shepard->smoothing_lengths[j];
    ++k;
  }
}

static void shepard_compute(void* context, 
                            int i, 
                            point_t* x,
                            real_t* values, 
                            vector_t* gradients)
{
  shepard_t* shepard = context;
  int N = shepard->N;

  // Compute the kernels and their gradients at x.
  real_t W[N];
  vector_t grad_W[N];
  shape_function_kernel_compute(shepard->W, shepard->xj, shepard->hj, N, x, W, grad_W);

  // Compute the values of the Shepard function.
  real_t sum_Wi = 0.0;
  for (i = 0; i < N; ++i)
    sum_Wi += W[i];
  if (sum_Wi == 0.0)
    memset(values, 0, sizeof(real_t) * N);
  else
  {
    for (i = 0; i < N; ++i)
      values[i] = W[i] / sum_Wi;
  }

  // Compute the gradients if needed.
  if (gradients != NULL)
  {
    if (sum_Wi == 0.0)
      memset(gradients, 0, sizeof(vector_t) * N);
    else
    {
      vector_t sum_grad_Wi = {.x = 0.0, .y = 0.0, .z = 0.0};
      for (i = 0; i < N; ++i)
      {
        sum_grad_Wi.x += grad_W[i].x;
        sum_grad_Wi.y += grad_W[i].y;
        sum_grad_Wi.z += grad_W[i].z;
      }

      // Use the quotient rule!
      for (i = 0; i < N; ++i)
      {
        gradients[i].x = (sum_Wi * grad_W[i].x - sum_grad_Wi.x * W[i]) / (sum_grad_Wi.x * sum_grad_Wi.x);
        gradients[i].y = (sum_Wi * grad_W[i].y - sum_grad_Wi.y * W[i]) / (sum_grad_Wi.y * sum_grad_Wi.y);
        gradients[i].z = (sum_Wi * grad_W[i].z - sum_grad_Wi.z * W[i]) / (sum_grad_Wi.z * sum_grad_Wi.z);
      }
    }
  }
}

static void shepard_dtor(void* context)
{
  shepard_t* shepard = context;
  polymec_free(shepard->xj);
  polymec_free(shepard->hj);
  polymec_free(shepard);
}

shape_function_t* shepard_shape_function_new(shape_function_kernel_t* kernel,
                                             point_cloud_t* domain,
                                             stencil_t* neighborhoods,
                                             real_t* smoothing_lengths)
{
  shepard_t* shepard = polymec_malloc(sizeof(shepard_t));
  shepard->W = kernel;
  shepard->domain = domain;
  shepard->neighborhoods = neighborhoods;
  shepard->smoothing_lengths = smoothing_lengths;

  // Count up the maximum neighborhood size and allocate storage.
  int max_neighborhood_size = -1;
  for (int i = 0; i < shepard->domain->num_points; ++i)
    max_neighborhood_size = MAX(max_neighborhood_size, stencil_size(shepard->neighborhoods, i));
  shepard->xj = polymec_malloc(sizeof(point_t) * max_neighborhood_size);
  shepard->hj = polymec_malloc(sizeof(real_t) * max_neighborhood_size);

  // Make sure our ghost points are consistent a representation of ghost points.
  stencil_exchange(shepard->neighborhoods, shepard->domain->points, 3, 0, MPI_REAL_T);
  stencil_exchange(shepard->neighborhoods, shepard->smoothing_lengths, 1, 0, MPI_REAL_T);

  shape_function_vtable vtable = {.neighborhood_size = shepard_neighborhood_size,
                                  .get_neighborhood_points = shepard_get_neighborhood_points,
                                  .set_neighborhood = shepard_set_neighborhood,
                                  .compute = shepard_compute,
                                  .dtor = shepard_dtor};
  return shape_function_new("Shepard", shepard, vtable);
}

