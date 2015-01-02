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
// the given imbalance tolerance. 
// NOTE: partitioning of point clouds using stencils instead of neighbor 
// NOTE: pairings is NOT supported, since stencils generally represent 
// NOTE: asymmetric neighbor relations, which produce directed graphs, and 
// NOTE: our graph partitioning algorithms all operate on undirected graphs.
exchanger_t* partition_point_cloud_with_neighbors(point_cloud_t** points, 
                                                  neighbor_pairing_t** neighbors, 
                                                  MPI_Comm comm, 
                                                  int* weights, 
                                                  real_t imbalance_tol);

#endif
