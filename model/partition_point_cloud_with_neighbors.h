// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARTITION_POINT_CLOUD_WITH_NEIGHBORS_H
#define POLYMEC_PARTITION_POINT_CLOUD_WITH_NEIGHBORS_H

#include "core/point_cloud.h"
#include "core/exchanger.h"
#include "model/neighbor_pairing.h"

// Given a global point cloud and a neighbor pairing connecting its points on 
// rank 0, this function partitions the points onto the processes on the given 
// communicator, replacing the point cloud and the stencil on each process 
// with the newly partitioned data. Each point can have a weight alloted to 
// it that characterizes its workload, and the load will be balanced within
// the given imbalance tolerance. A migrator capable of distributing data 
// to the partitions is returned.
// NOTE: partitioning of point clouds using stencils instead of neighbor 
// NOTE: pairings is NOT supported, since stencils generally represent 
// NOTE: asymmetric neighbor relations, which produce directed graphs, and 
// NOTE: our graph partitioning algorithms all operate on undirected graphs.
migrator_t* partition_point_cloud_with_neighbors(point_cloud_t** points, 
                                                  neighbor_pairing_t** neighbors, 
                                                  MPI_Comm comm, 
                                                  int* weights, 
                                                  real_t imbalance_tol);

#endif
