// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/partitioning.h"
#include "geometry/partition_point_cloud.h"
#include "model/partition_point_cloud_with_neighbors.h"

#if POLYMEC_HAVE_MPI
static void neighbor_pairing_distribute(neighbor_pairing_t** neighbors,
                                        MPI_Comm comm,
                                        int64_t* global_partition,
                                        size_t num_indices)
{
  START_FUNCTION_TIMER();

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Make sure we're all here.
  MPI_Barrier(comm);

  neighbor_pairing_t* global_pairing = *neighbors;
  neighbor_pairing_t* local_pairing = NULL;
  uint64_t vtx_dist[nprocs+1];
  if (rank == 0)
  {
    // Take stock of how many indices we'll have per process.
    int num_indices_p[nprocs];
    memset(num_indices_p, 0, sizeof(int) * nprocs);
    for (int i = 0; i < num_indices; ++i)
      num_indices_p[global_partition[i]]++;

    // Construct the distribution of indices ("vertices") for the partitioning.
    vtx_dist[0] = 0;
    for (int p = 0; p < nprocs; ++p)
      vtx_dist[p+1] = vtx_dist[p] + num_indices_p[p];

    // Now we create representations of pairings for each process. Since the
    // ordering of nodes (i, j) within a pair is unpredictable, we construct
    // them all simultaneously.
    log_debug("neighbor_pairing_distribute: creating pair representations on rank 0.");
    int num_pairs[nprocs];
    int_array_t* pairs[nprocs];
    exchanger_proc_map_t* sends[nprocs];
    exchanger_proc_map_t* receives[nprocs];
    int ghost_indices[nprocs];
    memset(num_pairs, 0, sizeof(int) * nprocs);
    memset(pairs, 0, sizeof(int_array_t*) * nprocs);
    memset(sends, 0, sizeof(exchanger_proc_map_t*) * nprocs);
    memset(receives, 0, sizeof(exchanger_proc_map_t*) * nprocs);
    memset(ghost_indices, 0, sizeof(int) * nprocs);
    int pos = 0, i, j;
    while (neighbor_pairing_next(global_pairing, &pos, &i, &j))
    {
      // Find the processes containing i and j.
      ASSERT(i >= vtx_dist[0]);
      ASSERT(i < vtx_dist[nprocs]);
      int pi = 0, pj = 0;
      while (vtx_dist[pi+1] < i) ++pi;
      while (vtx_dist[pj+1] < j) ++pj;

      // Deposit pairs in each of pi's and pj's arrays.
      if (pairs[pi] == NULL)
        pairs[pi] = int_array_new();
      int_array_append(pairs[pi], i);
      int_array_append(pairs[pi], j);
      if (pairs[pj] == NULL)
        pairs[pj] = int_array_new();
      int_array_append(pairs[pj], j);
      int_array_append(pairs[pj], i);

      // Set exchanger entries.
      if (pi != pj)
      {
        // Processes i and j must each have send/receive entries for this guy.

        // pi sends i to pj.
        if (sends[pi] == NULL)
          sends[pi] = exchanger_proc_map_new();
        exchanger_proc_map_add_index(sends[pi], pj, i);

        // pj sends j to pi.
        if (sends[pj] == NULL)
          sends[pj] = exchanger_proc_map_new();
        exchanger_proc_map_add_index(sends[pj], pi, j);

        // pi gets j from pj.
        if (receives[pi] == NULL)
          receives[pi] = exchanger_proc_map_new();
        exchanger_proc_map_add_index(receives[pi], pj, ghost_indices[pi]++);

        // pj gets i from pi.
        if (receives[pj] == NULL)
          receives[pj] = exchanger_proc_map_new();
        exchanger_proc_map_add_index(receives[pj], pi, ghost_indices[pj]++);
      }
    }

    // Now we create the pairings and ship them to the various processes.
    log_debug("neighbor_pairing_distribute: Distributing to other ranks.");
    serializer_t* ser = neighbor_pairing_serializer();
    byte_array_t* bytes = byte_array_new();
    for (int p = 0; p < nprocs; ++p)
    {
      // Share the vertex distribution with process p.
      if (p > 0)
        MPI_Send(vtx_dist, nprocs+1, MPI_UINT64_T, p, p, comm);

      // Set up the exchanger for process p's pairing.
      exchanger_t* ex = exchanger_new_with_rank(comm, p);
      if (sends[p] != NULL)
        exchanger_set_sends(ex, sends[p]);
      if (receives[p] != NULL)
        exchanger_set_receives(ex, receives[p]);

      // NOTE: The neighbor pairings consume the array data we've generated,
      // NOTE: so we call *_array_release_data_and_free() below to give it to
      // NOTE: them.
      neighbor_pairing_t* p_pairing = neighbor_pairing_new(global_pairing->name, num_pairs[p],
                                                           pairs[p]->data, ex);
      int_array_release_data_and_free(pairs[p]);

      if (p == 0)
      {
        // This here's our local pairing.
        local_pairing = p_pairing;
      }
      else
      {
        // Send the pairing to process p.
        size_t offset = 0;
        serializer_write(ser, p_pairing, bytes, &offset);
        log_debug("neighbor_pairing_distribute: sending neighbors to rank %d.", p);
        MPI_Send(&bytes->size, 1, MPI_INT, p, p, comm);
        MPI_Send(bytes->data, (int)bytes->size, MPI_BYTE, p, p, comm);

        // Clean up.
        byte_array_clear(bytes);
        neighbor_pairing_free(p_pairing);
      }
    }

    ser = NULL;
    byte_array_free(bytes);
  }
  else
  {
    // Receive the vertex distribution of the incoming pairing.
    MPI_Status status;
    MPI_Recv(vtx_dist, nprocs+1, MPI_UINT64_T, 0, rank, comm, &status);

    // Receive the size of the incoming pairing.
    int pairing_size;
    log_debug("neighbor_pairing_distribute: receiving neighbors from rank 0.");
    MPI_Recv(&pairing_size, 1, MPI_INT, 0, rank, comm, &status);

    // Now receive the pairing.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, pairing_size);
    MPI_Recv(bytes->data, pairing_size, MPI_BYTE, 0, rank, comm, &status);
    serializer_t* ser = neighbor_pairing_serializer();
    size_t offset = 0;
    local_pairing = serializer_read(ser, bytes, &offset);

    byte_array_free(bytes);
    ser = NULL;
  }

  *neighbors = local_pairing;

  // Clean up.
  if (global_pairing != NULL)
    neighbor_pairing_free(global_pairing);
  STOP_FUNCTION_TIMER();
}
#endif

bool partition_point_cloud_with_neighbors(point_cloud_t** points,
                                          neighbor_pairing_t** neighbors,
                                          MPI_Comm comm,
                                          int* weights,
                                          real_t imbalance_tol,
                                          point_cloud_field_t** fields,
                                          size_t num_fields)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);

  START_FUNCTION_TIMER();

#if POLYMEC_HAVE_MPI
  ASSERT((*points == NULL) || ((*points)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return true;

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
  int64_t* global_partition = NULL;
  if (rank == 0)
    log_debug("partition_point_cloud_with_neighbors: partitioning graph on rank 0.");

  global_partition = partition_graph(global_graph, comm, weights, imbalance_tol, true);
  if (global_partition == NULL)
    return false;

  // Break the neighbor pairing into chunks and send them to the other processes.
  neighbor_pairing_distribute(neighbors, comm, global_partition,
                              (*points != NULL) ? (*points)->num_points : 0);

  // Distribute the point cloud.
  distribute_point_cloud(points, comm, global_partition, fields, num_fields);

  // Clean up.
  if (global_graph != NULL)
    adj_graph_free(global_graph);
  if (global_partition != NULL)
    polymec_free(global_partition);

  STOP_FUNCTION_TIMER();

  return true;
#else
  return true;
#endif
}

