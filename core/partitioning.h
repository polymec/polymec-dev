// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PARTITIONING_H
#define POLYMEC_PARTITIONING_H

#include "core/adj_graph.h"
#include "core/exchanger.h"

/// \addtogroup core core
///@{

/// Partitions a (serial) global graph, creating and returning a global 
/// partition vector. The partitioning is only performed on rank 0. 
/// If broadcast is true, the partition vector is broadcasted and returned on 
/// all ranks. Otherwise only rank 0 returns the partition vector,
/// and all other ranks return NULL.
/// \collective Collective on comm.
int64_t* partition_graph(adj_graph_t* global_graph, 
                         MPI_Comm comm,
                         int* weights,
                         real_t imbalance_tol, 
                         bool broadcast);

/// Partitions a (serial) list of points, creating and returning a global 
/// partition vector. If broadcast is true, the partition vector is broadcasted
/// and returned on all ranks. Otherwise only rank 0 returns the partition vector,
/// and all other ranks return NULL.
/// \collective Collective on comm.
int64_t* partition_points(point_t* points,
                          size_t num_points,
                          MPI_Comm comm,
                          int* weights,
                          real_t imbalance_tol,
                          bool broadcast);

/// Repartitions a local graph, creating and returning a local partition 
/// vector with destination ranks included for ghost vertices.
/// \collective Collective on local_graph's communicator.
int64_t* repartition_graph(adj_graph_t* local_graph, 
                           exchanger_t* local_graph_ex,
                           int num_ghost_vertices,
                           int* weights,
                           real_t imbalance_tol);

/// Repartitions a local list of points, creating and returning a local 
/// partition vector with destination ranks included for ghost points.
/// \collective Collective on comm.
int64_t* repartition_points(point_t* local_points, 
                            size_t num_local_points,
                            MPI_Comm comm,
                            int* weights,
                            real_t imbalance_tol);

/// \struct redistribution
/// This struct contains information needed to redistribute data from the 
/// local process to each of its neighbors.
typedef struct
{
  int_array_t* send_procs;       // processes we sent data to
  int_array_t** send_indices;    // indices of the data we send to each proc
  int_array_t* receive_procs;    // processes we receive data from
  int_array_t** receive_indices; // indices of the data we receive from each proc
} redistribution_t;

/// Creates redistribution data from a local partition vector on the given 
/// MPI communicator for each calling rank. 
/// \memberof redistribution
/// \collective Collective on comm.
redistribution_t* redistribution_from_partition(MPI_Comm comm, 
                                                int64_t* local_partition,
                                                size_t num_vertices);

/// Destroys the given redistribution struct and all of its fields.
/// \memberof redistribution
void redistribution_free(redistribution_t* redist);

///@}

#endif
