// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/kd_tree.h"
#include "model/neighbor_pairing.h"

neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h);
neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h)
{
  // Toss all the points into a kd-tree.
  kd_tree_t* tree = kd_tree_new(cloud->points, cloud->num_points);

  // Now do a neighbor search -- everything that falls within h of a point 
  // is a neighbor.
  int_array_t* pairs = int_array_new();
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
      }
    }
    int_array_free(neighbors);
  }
  kd_tree_free(tree);

  exchanger_t* ex = exchanger_new(MPI_COMM_SELF);
  neighbor_pairing_t* pairing = neighbor_pairing_new("simple pairing", 
                                                     pairs->size/2,
                                                     pairs->data,
                                                     ex);
  int_array_release_data_and_free(pairs);
  return pairing;
}


