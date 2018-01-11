// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "core/adj_graph.h"
#include "core/exchanger.h"

#if POLYMEC_HAVE_MPI
#include "ptscotch.h"
#include "core/timer.h"
#endif

// This helper partitions a (serial) global graph, creating and returning a 
// global partition vector. It should only be called on rank 0. This is not 
// declared static because it may be used by other parts of polymec, but is 
// not part of the public API.
int64_t* partition_graph(adj_graph_t* global_graph, 
                         MPI_Comm comm,
                         int* weights,
                         real_t imbalance_tol);
int64_t* partition_graph(adj_graph_t* global_graph, 
                         MPI_Comm comm,
                         int* weights,
                         real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT(adj_graph_comm(global_graph) == MPI_COMM_SELF);

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT(rank == 0);

  SCOTCH_Dgraph dist_graph;
  SCOTCH_Num *vtx_weights = NULL, *xadj = NULL, *adj = NULL;
  int num_global_vertices = 0;
  if (rank == 0)
  {
    num_global_vertices = adj_graph_num_vertices(global_graph);

    // Extract the adjacency information.
    xadj = polymec_malloc(sizeof(SCOTCH_Num) * (num_global_vertices+1));
    int* edge_offsets = adj_graph_edge_offsets(global_graph);
    for (int i = 0; i <= num_global_vertices; ++i)
      xadj[i] = (SCOTCH_Num)edge_offsets[i];
    SCOTCH_Num num_arcs = xadj[num_global_vertices];
    adj = polymec_malloc(sizeof(SCOTCH_Num) * num_arcs);
    int* edges = adj_graph_adjacency(global_graph);
    for (int i = 0; i < xadj[num_global_vertices]; ++i)
      adj[i] = (SCOTCH_Num)edges[i];

    // Build a graph on rank 0.
    SCOTCH_dgraphInit(&dist_graph, MPI_COMM_SELF);
    if (weights != NULL)
    {
      vtx_weights = polymec_malloc(sizeof(SCOTCH_Num) * num_global_vertices);
      for (int i = 0; i < num_global_vertices; ++i)
        vtx_weights[i] = (SCOTCH_Num)weights[i];
    }
    SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_global_vertices, (SCOTCH_Num)num_global_vertices,
        xadj, NULL, vtx_weights, NULL, num_arcs, num_arcs,
        adj, NULL, NULL);
  }

  // Generate the global partition vector by scattering the global graph.
  int64_t* global_partition = NULL;
  if (rank == 0)
  {
    global_partition = polymec_malloc(sizeof(int64_t) * num_global_vertices);
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
    int result = SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    if (result != 0)
      polymec_error("Partitioning strategy could not be constructed.");
    result = SCOTCH_dgraphPart(&dist_graph, nprocs, &strategy, global_partition);
    if (result != 0)
      polymec_error("Partitioning failed.");
    SCOTCH_dgraphExit(&dist_graph);
    SCOTCH_stratExit(&strategy);

    if (vtx_weights != NULL)
      polymec_free(vtx_weights);
    polymec_free(xadj);
    polymec_free(adj);
  }

  // Return the global partition vector.
  STOP_FUNCTION_TIMER();
  return global_partition;
#else
  int num_global_vertices = adj_graph_num_vertices(global_graph);
  int64_t* P = polymec_malloc(sizeof(int64_t)*num_global_vertices);
  memset(P, 0, sizeof(int64_t)*num_global_vertices);
  return P;
#endif
}

// This helper repartitions a local graph, creating and returning a 
// local partition vector with destination ranks included for ghost cells.
int64_t* repartition_graph(adj_graph_t* local_graph, 
                           int num_ghost_vertices,
                           int* weights,
                           real_t imbalance_tol,
                           exchanger_t* local_graph_ex);
int64_t* repartition_graph(adj_graph_t* local_graph, 
                           int num_ghost_vertices,
                           int* weights,
                           real_t imbalance_tol,
                           exchanger_t* local_graph_ex)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm comm = adj_graph_comm(local_graph);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  int num_vertices = adj_graph_num_vertices(local_graph);

  // Extract the adjacency information.
  SCOTCH_Num* xadj = polymec_malloc(sizeof(SCOTCH_Num) * (num_vertices+1));
  int* edge_offsets = adj_graph_edge_offsets(local_graph);
  for (int i = 0; i <= num_vertices; ++i)
    xadj[i] = (SCOTCH_Num)edge_offsets[i];
  SCOTCH_Num num_arcs = xadj[num_vertices];
  SCOTCH_Num* adj = polymec_malloc(sizeof(SCOTCH_Num) * num_arcs);
  int* edges = adj_graph_adjacency(local_graph);
  for (int i = 0; i < xadj[num_vertices]; ++i)
    adj[i] = (SCOTCH_Num)edges[i];

  // Replace the ghost entries in adj with global indices.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  {
    index_t* global_indices = polymec_malloc(sizeof(index_t) * (num_vertices + num_ghost_vertices));
    for (int i = 0; i < num_vertices; ++i)
      global_indices[i] = (index_t)(vtx_dist[rank] + i);
    exchanger_exchange(local_graph_ex, global_indices, 1, 0, MPI_UINT64_T);
    for (int i = 0; i < num_arcs; ++i)
      adj[i] = global_indices[adj[i]];
    polymec_free(global_indices);
  }

  // Build a distributed graph.
  SCOTCH_Num* vtx_weights = NULL;
  if (weights != NULL)
  {
    vtx_weights = polymec_malloc(sizeof(SCOTCH_Num) * num_vertices);
    for (int i = 0; i < num_vertices; ++i)
      vtx_weights[i] = (SCOTCH_Num)weights[i];
  }
  SCOTCH_Dgraph dist_graph;
  SCOTCH_dgraphInit(&dist_graph, comm);
  SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_vertices, (SCOTCH_Num)num_vertices,
                     xadj, NULL, vtx_weights, NULL, num_arcs, num_arcs,
                     adj, NULL, NULL);

  // Generate the local partition vector by mapping the distributed graph.
  int64_t* local_partition = polymec_malloc(sizeof(int64_t) * (num_vertices + num_ghost_vertices));
  {
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATBALANCE;
    SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    int result = SCOTCH_dgraphPart(&dist_graph, nprocs, &strategy, local_partition);
    if (result != 0)
      polymec_error("Repartitioning failed.");
    SCOTCH_dgraphExit(&dist_graph);
    SCOTCH_stratExit(&strategy);
    polymec_free(adj);
    polymec_free(xadj);
  }
  if (vtx_weights != NULL)
    polymec_free(vtx_weights);

  // At this point, the local partition vector contains only the destination 
  // ranks of the local vertices. Now we talk amongst ourselves to figure 
  // out the destination ranks of the ghost vertices. This operation should 
  // preserve the ordering of the ghost vertices, too.
  exchanger_exchange(local_graph_ex, local_partition, 1, 0, MPI_UINT64_T);

  // Return the local partition vector.
  STOP_FUNCTION_TIMER();
  return local_partition;
#else
  int num_global_vertices = adj_graph_num_vertices(local_graph);
  int64_t* P = polymec_malloc(sizeof(int64_t)*num_global_vertices);
  memset(P, 0, sizeof(int64_t)*num_global_vertices);
  return P;
#endif
}

