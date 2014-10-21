// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "model/partition_point_cloud_with_neighbors.h"

// Defined in core/partition_mesh.c, but not part of API.
extern int64_t* partition_graph(adj_graph_t* global_graph, 
                                MPI_Comm comm, 
                                int* weights, 
                                real_t imbalance_tol);

// Defined in core/partition_point_cloud.c, but not part of API. 
extern void point_cloud_distribute(point_cloud_t** cloud, 
                                   MPI_Comm comm,
                                   int64_t* global_partition);

exchanger_t* partition_point_cloud_with_stencil(point_cloud_t** points, 
                                                stencil_t** stencil, 
                                                MPI_Comm comm, 
                                                int* weights, 
                                                real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);

#if POLYMEC_HAVE_MPI
  ASSERT((*points == NULL) || ((*points)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  // If points on rank != 0 are not NULL, we delete them.
  point_cloud_t* cloud = *points;
  if ((rank != 0) && (cloud != NULL))
  {
    point_cloud_free(cloud);
    *points = cloud = NULL; 
  }

  // Generate a global adjacency graph for the point cloud.
  adj_graph_t* global_graph = (cloud != NULL) ? graph_from_point_cloud_and_stencil(cloud, *stencil)
                                              : NULL;

#ifndef NDEBUG
  // Make sure there are enough points to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT((*points)->num_points > nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Reconstruct the stencil on the other domains.
  // FIXME

  // Distribute the point cloud.
  point_cloud_distribute(points, comm, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (cloud != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Clean up.
  if (global_graph != NULL)
    adj_graph_free(global_graph);
  if (global_partition != NULL)
    polymec_free(global_partition);

  // Return the migrator.
  return distributor;
#else
  return exchanger_new(comm);
#endif
}

exchanger_t* partition_point_cloud_with_neighbors(point_cloud_t** points, 
                                                  neighbor_pairing_t** neighbors, 
                                                  MPI_Comm comm, 
                                                  int* weights, 
                                                  real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);

#if POLYMEC_HAVE_MPI
  ASSERT((*points == NULL) || ((*points)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  // If points on rank != 0 are not NULL, we delete them.
  point_cloud_t* cloud = *points;
  if ((rank != 0) && (cloud != NULL))
  {
    point_cloud_free(cloud);
    *points = cloud = NULL; 
  }

  // Generate a global adjacency graph for the point cloud.
  adj_graph_t* global_graph = (cloud != NULL) ? graph_from_point_cloud_and_neighbors(cloud, *neighbors)
                                              : NULL;

#ifndef NDEBUG
  // Make sure there are enough points to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT((*points)->num_points > nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Reconstruct the neighbor pairing on the other domains.
  // FIXME

  // Distribute the point cloud.
  point_cloud_distribute(points, comm, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (cloud != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Clean up.
  if (global_graph != NULL)
    adj_graph_free(global_graph);
  if (global_partition != NULL)
    polymec_free(global_partition);

  // Return the migrator.
  return distributor;
#else
  return exchanger_new(comm);
#endif
}

