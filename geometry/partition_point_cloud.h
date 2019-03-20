// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARTITION_POINT_CLOUD_H
#define POLYMEC_PARTITION_POINT_CLOUD_H

#include "geometry/point_cloud_field.h"

// These functions provide partitioning and load balancing capabilities. In
// each of these, the imbalance tolerance is value between 0 and 1, defined as
// imbalance_tol = (process_work_load - balanced_work_load) / balanced_work_load.

/// \addtogroup geometry geometry
///@{

/// This function partitions the given point cloud on rank 0 with the given pointwise
/// load weights, distributing it to parallel subdomains on the given communicator
/// in such a way as to balance the load. The cloud on rank 0 (as
/// well as any non-NULL cloud on rank != 0) is consumed. In each case, the
/// cloud is replaced with a partitioned cloud. Any fields given are also partitioned
/// similarly. On success, this function returns true. If the per-process workload cannot
/// be balanced to within the imbalance tolerance, this function fails with no effect
/// and returns false.
/// \relates point_cloud
bool partition_point_cloud(point_cloud_t** cloud,
                           MPI_Comm comm,
                           int* weights,
                           real_t imbalance_tol,
                           point_cloud_field_t** fields,
                           size_t num_fields);

/// This function repartitions the given point cloud with the given load weights,
/// alloting them to parallel domains to balance their load. If weights is NULL, the
/// points are assigned equal weights. The cloud is consumed and replaced with a
/// repartitioned cloud, as are any given fields. On success, the function returns true.
/// If the per-process workload cannot be balanced to within the imbalance tolerance,
/// the function fails with no effect and returns false.
/// \relates point_cloud
bool repartition_point_cloud(point_cloud_t** cloud,
                             int* weights,
                             real_t imbalance_tol,
                             point_cloud_field_t** fields,
                             size_t num_fields);

//------------------------------------------------------------------------
// While partition_point_cloud and repartition_point_cloud are all-in-one
// point cloud partitioners, the following functions allow one to mix-n-match
// the pieces of the underlying algorithms.
//------------------------------------------------------------------------

/// Given a global partition vector, distribute the point cloud to the given communicator,
/// and distribute the data in the given fields accordingly. The cloud and fields are
/// replaced with partitioned equivalents.
/// \relates point_cloud
void distribute_point_cloud(point_cloud_t** cloud,
                            MPI_Comm comm,
                            int64_t* global_partition,
                            point_cloud_field_t** fields,
                            size_t num_fields);

/// Given a local partition vector, redistributes the point cloud and the given fields.
/// \relates point_cloud
void redistribute_point_cloud(point_cloud_t** cloud,
                              int64_t* local_partition,
                              point_cloud_field_t** fields,
                              size_t num_fields);

///@}

#endif

