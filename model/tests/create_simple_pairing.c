// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h);
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


