// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARTITION_POINT_CLOUD_H
#define POLYMEC_PARTITION_POINT_CLOUD_H

#include "core/point_cloud.h"
#include "core/exchanger.h"

// These functions provide partitioning and load balancing capabilities. In 
// each of these, the imbalance tolerance is value between 0 and 1, defined as 
// imbalance_tol = (process_work_load - balanced_work_load) / balanced_work_load.

// This function partitions the given point cloud on rank 0 with the given pointwise
// load weights, distributing it to parallel subdomains on the given communicator 
// in such a way as to balance the load. If weights is NULL, the points are assigned 
// equal weights. It creates and returns a migrator object that can be used 
// to distribute data from rank 0 to the partitions. The cloud on rank 0 (as 
// well as any non-NULL cloud on rank != 0) is consumed. In each case, the 
// cloud is replaced with a partitioned cloud. If the per-process workload 
// cannot be balanced to within the imbalance tolerance, this function fails 
// with no effect and returns NULL.
migrator_t* partition_point_cloud(point_cloud_t** cloud, 
                                  MPI_Comm comm, 
                                  int* weights, 
                                  real_t imbalance_tol);

// This function repartitions the given point cloud with the given load
// weights, alloting them to parallel domains to balance their load.
// If weights is NULL, the points are assigned equal weights.
// It creates and returns an migrator object that can be used to migrate 
// data from the old partition to the new. The cloud is consumed and replaced
// with a repartitioned cloud. If the per-process workload cannot be balanced 
// to within the imbalance tolerance, the function fails with no effect and 
// returns NULL.
migrator_t* repartition_point_cloud(point_cloud_t** cloud, int* weights, real_t imbalance_tol);

// While partition_point_cloud and repartition_point_cloud are all-in-one 
// point cloud partitioners, the following functions allow one to mix-n-match 
// the pieces of the underlying algorithms.

// This function creates a newly-allocated global partition vector that can be used 
// to distribute a global point cloud on rank 0 to all processes on the given communicator, 
// according to the given weights and the specified imbalance tolerance. Point clouds 
// on non-zero ranks are ignored.
int64_t* partition_vector_from_point_cloud(point_cloud_t* global_cloud, 
                                           MPI_Comm comm, 
                                           int* weights, 
                                           real_t imbalance_tol);

// Given a global partition vector, distribute the mesh to the given communicator, and 
// return a distributor object that can be used to distribute its data. The mesh is 
// replaced with a partitioned mesh.
migrator_t* distribute_point_cloud(point_cloud_t** cloud, MPI_Comm comm, int64_t* global_partition);

#endif

