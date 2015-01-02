// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/kd_tree.h"
#include "model/point_weight_function.h"
#include "model/neighbor_pairing.h"

static real_t hat_W(void* context, vector_t* y)
{
  real_t h = *((real_t*)context); // Spatial extent of W.
  real_t D = vector_mag(y);
  return (D < h) ? 1.0 - D/h : 0.0;
}

static vector_t hat_grad(void* context, vector_t* y)
{
  real_t h = *((real_t*)context); // Spatial extent of W.
  real_t D = vector_mag(y);
  static vector_t zero = {.x = 0.0, .y = 0.0, .z = 0.0};
  vector_t g = {.x = -y->x/(D*h), .y = -y->y/(D*h), .z = -y->z/(D*h)};
  return (D < h) ? g : zero;
}

static void hat_dtor(void* context)
{
  polymec_free(context);
}

static point_weight_function_t* hat_function_new(real_t h)
{
  real_t* Wh = polymec_malloc(sizeof(real_t));
  *Wh = h;
  point_weight_function_vtable vtable = {.value = hat_W, 
                                         .gradient = hat_grad, 
                                         .dtor = hat_dtor};
  return point_weight_function_new("Hat function", Wh, vtable);
}

neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h)
{
  point_weight_function_t* W = hat_function_new(h);

  // Toss all the points into a kd-tree.
  kd_tree_t* tree = kd_tree_new(cloud->points, cloud->num_points);

  // Now do a neighbor search -- everything that falls within h of a point 
  // is a neighbor.
  int_array_t* pairs = int_array_new();
  real_array_t* weights = real_array_new();
  for (int i = 0; i < cloud->num_points; ++i)
  {
    point_t* xi = &cloud->points[i];
    int_array_t* neighbors = kd_tree_within_radius(tree, &cloud->points[i], h);
    for (int j = 0; j < neighbors->size; ++j)
    {
      int k = neighbors->data[j];
      if (k > i)
      {
        int_array_append(pairs, i);
        int_array_append(pairs, k);
        point_t* xk = &cloud->points[k];
        vector_t y;
        point_displacement(xk, xi, &y);
        real_array_append(weights, point_weight_function_value(W, &y));
      }
    }
    int_array_free(neighbors);
  }
  kd_tree_free(tree);
  point_weight_function_free(W);

  exchanger_t* ex = exchanger_new(MPI_COMM_SELF);
  neighbor_pairing_t* pairing = neighbor_pairing_new("simple pairing", 
                                                     pairs->size/2,
                                                     pairs->data,
                                                     weights->data,
                                                     ex);
  int_array_release_data_and_free(pairs);
  real_array_release_data_and_free(weights);
  return pairing;
}


