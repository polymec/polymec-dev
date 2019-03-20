// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARTITION_MESH_H
#define POLYMEC_PARTITION_MESH_H

#include "geometry/polymesh_field.h"

/// \addtogroup geometry geometry
///@{

/// This function partitions the given polymesh on rank 0 with the given cell-centered
/// load weights, distributing the cells to parallel domains on the given
/// communicator to balance the load. If weights is NULL, the cells are assumed
/// all to have equal weights. The mesh is replaced with a partitioned mesh.
/// Any fields given are also partitioned similarly. On success, this function
/// returns true. If the per-process workload cannot be balanced to within the
/// imbalance tolerance, this function fails with no effect and returns false.
/// \relates polymesh
bool partition_polymesh(polymesh_t** mesh,
                        MPI_Comm comm,
                        int* weights,
                        real_t imbalance_tol,
                        polymesh_field_t** fields,
                        size_t num_fields);

/// This function repartitions the given mesh with the given load weights,
/// alloting the cells to parallel domains to balance their load. If weights is
/// NULL, the cells are all assumed to have equal weights. The mesh is replaced
/// with a partitioned mesh. Any fields given are also partitioned similarly.
/// On success, this function returns true. If the per-process workload cannot
/// be balanced to within the imbalance tolerance, this function fails with no
/// effect and returns false.
/// \relates polymesh
bool repartition_polymesh(polymesh_t** mesh,
                          int* weights,
                          real_t imbalance_tol,
                          polymesh_field_t** fields,
                          size_t num_fields);

// While partition_polymesh and repartition_polymesh are all-in-one mesh partitioners, the
// following functions allow one to mix-n-match the pieces of the underlying algorithms.

/// This function creates a newly-allocated global partition vector that can be used
/// to distribute a global polymesh on rank 0 to all processes on the given communicator,
/// according to the given weights and the specified imbalance tolerance. Mesh objects
/// on non-zero ranks are ignored. If broadcast is set to true, the global
/// partition vector is returned on every process--otherwise, it is just
/// returned on rank 0 of the given communicator.
/// \relates polymesh
int64_t* partition_vector_from_polymesh(polymesh_t* global_mesh,
                                        MPI_Comm comm,
                                        int* weights,
                                        real_t imbalance_tol,
                                        bool broadcast);

/// Given a global partition vector, distributes the mesh from rank 0 to all
/// processes in the given communicator according to the global partition vector,
/// distributing field data accordingly. The mesh and fields are replaced with
/// partitioned equivalents.
/// \relates polymesh
void distribute_polymesh(polymesh_t** mesh,
                         MPI_Comm comm,
                         int64_t* global_partition,
                         polymesh_field_t** fields,
                         size_t num_fields);

/// Given a local partition vector, redistributes the mesh and the given fields.
/// \relates polymesh
void redistribute_polymesh(polymesh_t** mesh,
                           int64_t* local_partition,
                           polymesh_field_t** fields,
                           size_t num_fields);

///@}

#endif

