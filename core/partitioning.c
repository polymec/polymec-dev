// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "core/partitioning.h"
#include "core/hilbert.h"
#include "core/parallel_sort.h"
#include "core/timer.h"

#if POLYMEC_HAVE_MPI
#include "ptscotch.h"
#endif

#if POLYMEC_HAVE_MPI

// This helper creates a Hilbert space-filling curve from a set of points.
static hilbert_t* hilbert_from_points(point_t* points, size_t num_points)
{
  bbox_t bbox = {.x1 = REAL_MAX, .x2 = -REAL_MAX, 
                 .y1 = REAL_MAX, .y2 = -REAL_MAX, 
                 .z1 = REAL_MAX, .z2 = -REAL_MAX};
  for (int i = 0; i < num_points; ++i)
    bbox_grow(&bbox, &points[i]);

  // (Handle lower dimensional point distributions gracefully.)
  if (reals_equal(bbox.x2, bbox.x1))
  {
    bbox.x1 -= 0.5;
    bbox.x2 += 0.5;
  }
  if (reals_equal(bbox.y2, bbox.y1))
  {
    bbox.y1 -= 0.5;
    bbox.y2 += 0.5;
  }
  if (reals_equal(bbox.z2, bbox.z1))
  {
    bbox.z1 -= 0.5;
    bbox.z2 += 0.5;
  }
  return hilbert_new(&bbox);
}

static void extract_graph_info(adj_graph_t* graph,
                               int* weights,
                               SCOTCH_Num** xadj, 
                               SCOTCH_Num** adj, 
                               SCOTCH_Num *num_arcs,
                               SCOTCH_Num** vtx_weights)
{
  size_t num_vertices = adj_graph_num_vertices(graph);
  *xadj = polymec_malloc(sizeof(SCOTCH_Num) * (num_vertices+1));
  int* edge_offsets = adj_graph_edge_offsets(graph);
  for (int i = 0; i <= num_vertices; ++i)
    (*xadj)[i] = (SCOTCH_Num)edge_offsets[i];
  *num_arcs = (*xadj)[num_vertices];
  *adj = polymec_malloc(sizeof(SCOTCH_Num) * (*num_arcs));
  int* edges = adj_graph_adjacency(graph);
  for (int i = 0; i < *num_arcs; ++i)
    (*adj)[i] = (SCOTCH_Num)edges[i];
  if (weights != NULL)
  {
    *vtx_weights = polymec_malloc(sizeof(SCOTCH_Num) * num_vertices);
    for (int i = 0; i < num_vertices; ++i)
      (*vtx_weights)[i] = (SCOTCH_Num)weights[i];
  }
}
#endif

int64_t* partition_graph(adj_graph_t* global_graph, 
                         MPI_Comm comm,
                         int* weights,
                         real_t imbalance_tol,
                         bool broadcast)
{
#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");
  START_FUNCTION_TIMER();

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((global_graph != NULL) || (rank != 0));
  ASSERT((global_graph == NULL) || (adj_graph_comm(global_graph) == MPI_COMM_SELF));
  ASSERT((rank == 0) || broadcast);

  SCOTCH_Dgraph dist_graph;
  SCOTCH_Num *vtx_weights = NULL, *xadj = NULL, *adj = NULL;
  size_t num_global_vertices = 0;
  if (rank == 0)
  {
    num_global_vertices = adj_graph_num_vertices(global_graph);

    // Extract the adjacency information.
    SCOTCH_Num num_arcs;
    extract_graph_info(global_graph, weights, &xadj, &adj, &num_arcs, &vtx_weights);

    // Build a graph on rank 0.
    SCOTCH_dgraphInit(&dist_graph, MPI_COMM_SELF);
    SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_global_vertices, (SCOTCH_Num)num_global_vertices,
                       xadj, NULL, vtx_weights, NULL, num_arcs, num_arcs, adj, NULL, NULL);
#ifndef NDEBUG
    SCOTCH_dgraphCheck(&dist_graph);
#endif
  }

  // Cut up the graph and generate the global partition vector.
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

  // Now broadcast the partition vector if we're asked to.
  if (broadcast)
  {
    MPI_Bcast(&num_global_vertices, 1, MPI_SIZE_T, 0, comm);
    if (rank != 0)
      global_partition = polymec_malloc(sizeof(int64_t) * num_global_vertices);
    MPI_Bcast(global_partition, (int)num_global_vertices, MPI_INT64_T, 0, comm);
  }

  // Return the global partition vector.
  STOP_FUNCTION_TIMER();
  return global_partition;
#else
  size_t num_global_vertices = adj_graph_num_vertices(global_graph);
  int64_t* P = polymec_calloc(num_global_vertices, sizeof(int64_t));
  return P;
#endif
}

int64_t* partition_graph_n_ways(adj_graph_t* global_graph, 
                                size_t n,
                                int* weights,
                                real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");
  ASSERT(global_graph != NULL);

  START_FUNCTION_TIMER();
  size_t num_global_vertices = adj_graph_num_vertices(global_graph);
  int64_t* global_partition = polymec_calloc(num_global_vertices, sizeof(int64_t));

  if (n == 1)
  {
    // We're finished!
    STOP_FUNCTION_TIMER();
    return global_partition;
  }

  // Extract the adjacency information.
  SCOTCH_Num *vtx_weights = NULL, *xadj = NULL, *adj = NULL, num_arcs = 0;
  extract_graph_info(global_graph, weights, &xadj, &adj, &num_arcs, &vtx_weights);

  // Build a graph.
  SCOTCH_Dgraph dist_graph;
  int result = SCOTCH_dgraphInit(&dist_graph, MPI_COMM_SELF);
  if (result != 0)
    polymec_error("partition_graph_n_ways: Graph could not be initialized.");
  result = SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_global_vertices, 
                              (SCOTCH_Num)num_global_vertices, xadj, NULL, 
                              vtx_weights, NULL, num_arcs, num_arcs, adj, NULL, NULL);
  if (result != 0)
    polymec_error("partition_graph_n_ways: Graph could not be constructed.");
#ifndef NDEBUG
  SCOTCH_dgraphCheck(&dist_graph);
#endif

  // Cut up the graph -> global partition vector.
  SCOTCH_Strat strategy;
  SCOTCH_stratInit(&strategy);
  SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
  result = SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, 1, n, (double)imbalance_tol);
  if (result != 0)
    polymec_error("Partitioning strategy could not be constructed.");
  log_debug("partition_graph_n_ways: Partitioning strategy constructed for %d parts.", n);
  result = SCOTCH_dgraphPart(&dist_graph, n, &strategy, global_partition);
  if (result != 0)
    polymec_error("Partitioning failed.");
  log_debug("partition_graph_n_ways: Computed partition vector.");
  SCOTCH_dgraphExit(&dist_graph);
  SCOTCH_stratExit(&strategy);
  if (vtx_weights != NULL)
    polymec_free(vtx_weights);
  polymec_free(xadj);
  polymec_free(adj);
  // Return the global partition vector.
  STOP_FUNCTION_TIMER();
  return global_partition;
#else
  // Do a naive partitioning.
  log_urgent("partition_graph_n_ways: MPI disabled. Performing naive partitioning. This probably isn't what you want.");
  size_t num_global_vertices = adj_graph_num_vertices(global_graph);
  int64_t* P = polymec_malloc(sizeof(int64_t)*num_global_vertices);
  int n_per_piece = (int)(num_global_vertices / n);
  for (int i = 0; i < num_global_vertices; ++i)
    P[i] = (int64_t)(i / n_per_piece);
  return P;
#endif
}

int64_t* partition_points(point_t* points,
                          size_t num_points,
                          MPI_Comm comm,
                          int* weights,
                          real_t imbalance_tol,
                          bool broadcast)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (points != NULL));

  int64_t* global_partition = NULL;

  // Do an n-way partitioning on rank 0.
  if (rank == 0)
  {
    global_partition = partition_points_n_ways(points, num_points, nprocs,
                                               weights, imbalance_tol);
  }

  // Now broadcast the partition vector if we're asked to.
  if (broadcast)
  {
    MPI_Bcast(&num_points, 1, MPI_SIZE_T, 0, comm);
    if (rank != 0)
      global_partition = polymec_malloc(sizeof(int64_t) * num_points);
    MPI_Bcast(global_partition, (int)num_points, MPI_INT64_T, 0, comm);
  }

  STOP_FUNCTION_TIMER();
  return global_partition;

#else
  // This is dumb, but we were asked for it.
  int64_t* global_partition = polymec_calloc(num_points, sizeof(int64_t));
  return global_partition;
#endif
}

// This helper is a comparison function used to sort (Hilbert index, weight) tuples.
// Only the Hilbert index factors into the ordering.
static int hilbert_comp(const void* l, const void* r)
{
  const index_t* li = l;
  const index_t* ri = r;
  return (li[1] < ri[1]) ? -1
                         : (li[1] > ri[1]) ? 1
                                           : 0;
}

int64_t* partition_points_n_ways(point_t* points,
                                 size_t num_points,
                                 size_t n,
                                 int* weights,
                                 real_t imbalance_tol)
{
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));

  // Handle the trivial case.
  if (n == 1)
  {
    int64_t* global_partition = polymec_calloc(num_points, sizeof(int64_t));
    STOP_FUNCTION_TIMER();
    return global_partition;
  }

#ifndef NDEBUG
  // Make sure there are enough points to go around for the number we're given.
  ASSERT(num_points > n);
#endif

  // Set up a Hilbert space filling curve that can map the given points to indices.
  bbox_t bbox;
  bbox_make_empty_set(&bbox);
  for (size_t i = 0; i < num_points; ++i)
    bbox_grow(&bbox, &points[i]);
  hilbert_t* hilbert = hilbert_new(&bbox);

  // Create an array of 3-tuples containing the 
  // (point index, Hilbert index, weight) of each point. 
  // Partitioning the points amounts to sorting this array and breaking it into parts whose 
  // work is equal. Also sum up the work on the points.
  index_t* part_array = polymec_malloc(sizeof(index_t) * 3 * num_points);
  uint64_t total_work = 0;
  if (weights != NULL)
  {
    for (int i = 0; i < num_points; ++i)
    {
      part_array[3*i]   = i;
      part_array[3*i+1] = hilbert_index(hilbert, &points[i]);
      part_array[3*i+2] = (index_t)weights[i];
      total_work += weights[i];
    }
  }
  else
  {
    for (int i = 0; i < num_points; ++i)
    {
      part_array[3*i]   = i;
      part_array[3*i+1] = hilbert_index(hilbert, &points[i]);
      part_array[3*i+2] = 1;
    }
    total_work = num_points;
  }

  // Sort the array.
  qsort(part_array, num_points, 3*sizeof(index_t), hilbert_comp);
  release_ref(hilbert);

  // Now we need to break it into parts of equal work.
  real_t work_per_proc = 1.0 * total_work / n;
  int part_offsets[n+1];
  real_t part_work[n+1];
  part_offsets[0] = 0;
  part_work[0] = 0.0;
  for (int p = 0; p < n; ++p)
  {
    int i = part_offsets[p];
    real_t work = 0.0, last_weight = 0.0, cum_work = part_work[p];
    while ((cum_work < ((p+1) * work_per_proc)) && (i < num_points))
    {
      last_weight = 1.0 * part_array[3*i+2];
      work += last_weight;
      cum_work += last_weight;
      ++i;
    }

    // If we've obviously overloaded this process, back up one step.
    if (((work - work_per_proc)/work_per_proc > imbalance_tol) && 
        ((work_per_proc - (work - last_weight) <= imbalance_tol)))
    {
      --i;
      work -= last_weight;
    }
    part_offsets[p+1] = i;
    part_work[p+1] = cum_work;
  }
   
  // Now we create the global partition vector and fill it.
  int64_t* global_partition = polymec_malloc(sizeof(int64_t) * num_points);
  int k = 0;
  for (int p = 0; p < n; ++p)
  {
    for (int i = part_offsets[p]; i < part_offsets[p+1]; ++i, ++k)
    {
      index_t j = part_array[3*k];
      global_partition[j] = p;
    }
  }
    
  // Clean up.
  polymec_free(part_array);

  STOP_FUNCTION_TIMER();
  return global_partition;
}

int64_t* repartition_graph(adj_graph_t* local_graph, 
                           exchanger_t* local_graph_ex,
                           int num_ghost_vertices,
                           int* weights,
                           real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm comm = adj_graph_comm(local_graph);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Extract the adjacency information.
  SCOTCH_Num *vtx_weights = NULL, *xadj = NULL, *adj = NULL, num_arcs = 0;
  extract_graph_info(local_graph, weights, &xadj, &adj, &num_arcs, &vtx_weights);

  // Replace the ghost entries in adj with global indices.
  size_t num_vertices = adj_graph_num_vertices(local_graph);
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
  size_t num_global_vertices = adj_graph_num_vertices(local_graph);
  int64_t* P = polymec_calloc(num_global_vertices, sizeof(int64_t));
  return P;
#endif
}

#if POLYMEC_HAVE_MPI 

// This rebalances the workload on the given set of points, changing the distribution 
// of the points on the processors but preserving their ordering. The points are 
// assumed to be 4-wide sets of indices as ordered within repartition_points
// below. This procedure is serial, so the function must be used sparingly. 
// It returns true if it can balance the loads to within the imbalance toleranace,
// false if not.
static bool balance_loads(MPI_Comm comm, 
                          real_t imbalance_tol,
                          index_t** local_array,
                          size_t* local_array_size)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  index_t* array = *local_array;
  size_t array_size = *local_array_size;

  // Sum the total load on all processes, and find the ideal load per process.
  real_t my_load = 0.0;
  for (int i = 0; i < array_size; ++i)
    my_load += 1.0 * array[4*i+3];
  log_debug("balance_loads: Current process workload is %g.", my_load);

  real_t total_load;
  MPI_Allreduce(&my_load, &total_load, 1, MPI_REAL_T, MPI_SUM, comm);
  real_t ideal_proc_load = total_load / nprocs;
  log_debug("balance_loads: Ideal process workload is %g.", ideal_proc_load);

  // Are we balanced already?
  {
    real_t max_proc_load, min_proc_load;
    MPI_Allreduce(&my_load, &max_proc_load, 1, MPI_REAL_T, MPI_MAX, comm);
    MPI_Allreduce(&my_load, &min_proc_load, 1, MPI_REAL_T, MPI_MIN, comm);
    real_t imbalance = MAX(ABS(max_proc_load - total_load/nprocs), 
                           ABS(min_proc_load - total_load/nprocs));
    log_debug("balance_loads: Current maximum imbalance is %g%%.", imbalance * 100.0 * nprocs / total_load);
    if (imbalance <= imbalance_tol * total_load / nprocs)
      return true;
  }

  // Now loop through the processes and carve off the ideal workload.
  int balanced = 1; // let's be optimistic!
  for (int p = 0; p < nprocs; ++p)
  {
    if (rank == p)
    {
      if (rank > 0)
      {
        // Catch the number of points we'll be getting from process (rank - 1).
        int num_points;
        MPI_Recv(&num_points, 1, MPI_INT, rank-1, 0, comm, MPI_STATUS_IGNORE);

        if (num_points > 0)
        {
          log_debug("balance_loads: Receiving %d points from previous process.", num_points);
          // Now catch the data for each of these points.
          index_t buffer[4*num_points];
          MPI_Recv(buffer, 4*num_points, MPI_INDEX_T, rank-1, 0, comm, MPI_STATUS_IGNORE);

          // Prepend these indices/loads to our local array.
          array_size += num_points;
          *local_array_size = array_size;
          *local_array = polymec_realloc(*local_array, 4 * sizeof(index_t) * array_size);
          array = *local_array;
          memmove(&(array[4*num_points]), array, 4*sizeof(index_t) * (array_size - num_points));
          memcpy(array, buffer, 4*sizeof(index_t) * num_points);
        }
        else if (num_points < 0)
        {
          // If the number of points we receive is negative, the datum is 
          // -work needed by process (rank - 1). So we count off (from the 
          // beginning of our array) the number of points that most closely 
          // satisfies this quantity, and send back that number of points.
          int work_needed = -num_points; 
          int num_points_sent = 0, work_given = 0, delta = 0;
          while (work_given <= work_needed)
          {
            delta = (int)array[4*num_points_sent+3];
            work_given += delta;
            ++num_points_sent;
          }

          // Step back a point if it makes sense.
          if ((work_given > work_needed) && 
              (abs(work_given - delta - work_needed) < (work_given - work_needed)))
          {
            work_given -= delta;
            --num_points_sent;
          }

          // Send the number of points back.
          log_debug("balance_loads: Sending %d points to previous process.", num_points_sent);
          MPI_Send(&num_points_sent, 1, MPI_INT, rank-1, 0, comm);

          // Send the data for those points.
          MPI_Send(array, 4*num_points_sent, MPI_INDEX_T, rank-1, 0, comm);

          // If we sent back more work than was strictly requested, wait to 
          // see whether it was rejected. A 1 response means accepted, 0 
          // rejected. If rejected, we proceed having passed one fewer 
          // point back.
          int accepted;
          MPI_Recv(&accepted, 1, MPI_INT, rank-1, 0, comm, MPI_STATUS_IGNORE);
          if (!accepted)
            --num_points_sent;

          // Remove the points we've sent from the front of our array.
          if (num_points_sent > 0)
          {
            memmove(array, &array[4*num_points_sent], 4*sizeof(index_t) * num_points_sent);
            array_size -= num_points_sent;
            *local_array_size = array_size;
          }
        }
      }

      // March through our points and add up the workload.
      real_t load = 0, delta;
      size_t num_points_kept = 0;
      while ((load < ideal_proc_load) && (num_points_kept < array_size))
      {
        delta = (real_t)array[4*num_points_kept+3];
        load += delta;
        ++num_points_kept;
      }

      // If there's a process "to our right," dump any excess work there.
      real_t excess_load = load - ideal_proc_load;
      real_t imbalance = excess_load / ideal_proc_load;
      if ((excess_load >= 0.0) || (num_points_kept < array_size))
      {
        if (rank == (nprocs - 1)) // last process!
        {
          // Can we fit the extra work on the process and still fall under 
          // our load imbalance tolerance?
          if (imbalance > imbalance_tol)
          {
            log_debug("balance_loads: Load on last process (%d) is unbalanced (%g%%).", 
                      (int)load, imbalance * 100.0);
            balanced = 0; 
            break;
          }
          else
          {
            // Just keep all the points on this process.
            if (excess_load > 0)
              log_debug("balance_loads: Keeping excess load on last process.");
            num_points_kept = array_size;
          }
        }
        else // We have a sucker to pass our extra work to!
        {
          size_t num_points_sent = array_size - num_points_kept;

          // Tell the sucker how many points we're giving him.
          log_debug("balance_loads: Sending %d points to next process.", num_points_sent);
          MPI_Send(&num_points_sent, 1, MPI_INT, rank+1, 0, comm);

          if (num_points_sent > 0)
          {
            // Pass the data along.
            MPI_Send(&array[4*num_points_kept], (int)(4*num_points_sent), MPI_INDEX_T, 
                     rank+1, 0, comm);
            array_size -= num_points_sent;
            *local_array = polymec_realloc(*local_array, 4 * sizeof(index_t) * array_size);
            *local_array_size = array_size;
          }
          else
          {
            // Just keep all the points on this process.
            num_points_kept = array_size;
          }
        }
      }
      else if (excess_load < 0.0) 
      {
        // Can we grab some work from our right-ward neighbor?
        if (rank < (nprocs - 1)) 
        {
          // We send the negation of how much work we would like to balance our 
          // load.
          int needed_work = (int)(ceil(ideal_proc_load - load));
          needed_work *= -1; // Flip the sign to signify that we want work.
          MPI_Send(&needed_work, 1, MPI_INT, rank+1, 0, comm);

          // Get back the number of points we're receiving.
          int num_points_received;
          MPI_Recv(&num_points_received, 1, MPI_INT, rank+1, 0, comm, MPI_STATUS_IGNORE);
          log_debug("balance_loads: Receiving %d points from next process.", num_points_received);

          // Get data for these points.
          index_t buffer[4*num_points_received];
          MPI_Recv(buffer, 4*num_points_received, MPI_INDEX_T, rank+1, 0, comm, MPI_STATUS_IGNORE);

          // Are we now overloaded? If so, return a 0 to signify that we 
          // reject the last of the points given, and a 1 to signify that 
          // we accept the lot.
          int extra_work = 0, response = 1;
          for (int i = 0; i < num_points_received; ++i)
            extra_work += (int)buffer[4*i+3];
#ifndef NDEBUG
          if (num_points_received > 0)
          {
            int last_load = (int)buffer[4*(num_points_received-1)+3];
            ASSERT(extra_work - last_load + load <= ideal_proc_load); // safety guarantee!
          }
#endif
          if (extra_work + load > ideal_proc_load)
          {
            real_t new_imbalance = (1.0 * (load + excess_load) - ideal_proc_load) / ideal_proc_load;
            if (new_imbalance > imbalance_tol)
            {
              response = 0; // we reject the last point.
              --num_points_received;
              log_debug("balance_loads: Rejecting final point from next process.");
            }
          }

          // Send the response.
          MPI_Send(&response, 1, MPI_INT, rank+1, 0, comm);

          // Now append the extra points to our local array.
          array_size += num_points_received;
          *local_array_size = array_size;
          *local_array = polymec_realloc(*local_array, 4 * sizeof(index_t) * array_size);
          array = *local_array;
          memcpy(&array[4*(array_size-num_points_received)], buffer, 
                 4 * sizeof(index_t) * num_points_received);
        }
        else if (ABS(imbalance) > imbalance_tol) // we are the last process
        {
          log_debug("balance_loads: Load on last process (%d) is unbalanced (%g%%).", 
                    (int)load, imbalance * 100.0);
          balanced = 0; 
          break;
        }
      }
    }
  }

  // Let's all agree on whether the load was balanced.
  MPI_Allreduce(MPI_IN_PLACE, &balanced, 1, MPI_INT, MPI_MIN, comm);
  if (balanced == 1)
  {
    log_debug("balance_loads: Targeting %d local points.", array_size);
    log_debug("balance_loads: Successfully balanced workloads to within %d%%.", 
              (int)(round(imbalance_tol * 100.0)));
  }
  else
  {
    log_debug("balance_loads: Could not balance workloads to within %d%%.",
              (int)(round(imbalance_tol * 100.0)));
  }

  STOP_FUNCTION_TIMER();
  return (balanced == 1);
}

// This creates local partition and load vectors using the information in the sorted 
// distributed array. This is to be used only by repartition_points().
static int64_t* partition_vector_from_balanced_array(MPI_Comm comm, 
                                                     index_t* balanced_array, 
                                                     size_t balanced_array_size,
                                                     size_t orig_array_size)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Find out who is sending us data by reading off the ranks of the processes that 
  // originally owned each point.

  // Store the number of points we are receiving from rank p in 
  // num_points_from_rank[p] and those we are sending to rank p 
  // in num_points_to_rank[p].
  int num_points_from_rank[nprocs], num_points_to_rank[nprocs];
  {
    memset(num_points_from_rank, 0, sizeof(int) * nprocs);
    for (int i = 0; i < balanced_array_size; ++i)
    {
      // We read off the rank that originally owned this point.
      int point_rank = (int)balanced_array[4*i+1];

      // If it's not this rank, we are expecting it from its owner.
      if (rank != point_rank)
        ++num_points_from_rank[point_rank];
    }

    // Convey to everyone how many points they are sending to whom.
    MPI_Alltoall(num_points_from_rank, 1, MPI_INT, num_points_to_rank, 1, MPI_INT, comm);

    if (log_level() == LOG_DEBUG) // Throw the poor dev a bone.
    {
      for (int p = 0; p < nprocs; ++p)
      {
        if (num_points_from_rank[p] > 0) 
          log_debug("partition_vector_from_balanced_array: Getting %d points from rank %d.", num_points_from_rank[p], p);
        if (num_points_to_rank[p] > 0)
          log_debug("partition_vector_from_balanced_array: Assigning %d points to rank %d.", num_points_to_rank[p], p);
      }
    }
  }

  // Post receives/sends for all nonzero sets of points. Of course, since we 
  // are given only the information about who we're receiving points from 
  // (during the migration process, that is), we need to have the receiving 
  // processes send the sending processes those points that they are supposed 
  // to be, uhm, sending. So this might seem backwards.

  // Do a preliminary count of processes that interact with this one. Note
  // that "sends" and "receives" here refer to the MPI_Sends and MPI_Receives
  // below, which run in the opposite direction of the pending migration.
  // ("Hmmmmm....")
  int num_sends = 0, num_receives = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    if (num_points_from_rank[p] > 0)
      ++num_sends;

    if (num_points_to_rank[p] > 0)
      ++num_receives;
  }

  int_ptr_unordered_map_t* sends = int_ptr_unordered_map_new();
  if (num_sends + num_receives > 0)
  {
    // Now send the indices of the points that the receivers need in order to construct
    // a partition vector.
    int i_req = 0;
    MPI_Request requests[num_receives + num_sends];
    int_ptr_unordered_map_t* receives = int_ptr_unordered_map_new();
    for (int p = 0; p < nprocs; ++p)
    {
      // Post receives on the sending process end.
      if (num_points_to_rank[p] > 0) 
      {
        int* indices = polymec_malloc(sizeof(int) * num_points_to_rank[p]);
        int_ptr_unordered_map_insert_with_v_dtor(sends, p, indices, polymec_free);
        MPI_Irecv(indices, num_points_to_rank[p], MPI_INT, p, 0, comm, &requests[i_req]);
        ++i_req;
      }

      // Start sends on the receiving process end.
      if (num_points_from_rank[p] > 0) 
      {
        int* indices = polymec_malloc(sizeof(int) * num_points_from_rank[p]);
        int j = 0;
        for (int i = 0; i < balanced_array_size; ++i)
        {
          int rank1 = (int)balanced_array[4*i+1];
          if (rank1 == p)
            indices[j++] = (int)balanced_array[4*i+2];
        }
        int_ptr_unordered_map_insert_with_v_dtor(receives, p, indices, polymec_free);
        MPI_Isend(indices, num_points_from_rank[p], MPI_INT, p, 0, comm, &requests[i_req]);
        ++i_req;
      }
    }
    ASSERT(i_req == num_sends + num_receives);

    // Wait for the messages to be passed.
    MPI_Status statuses[num_receives + num_sends];
    MPI_Waitall(num_receives + num_sends, requests, statuses);

    int_ptr_unordered_map_free(receives);
  }

  // At this point, we have the indices we need to send to each process.
  // We can now construct a partition vector whose values for these indices
  // will be the corresponding process rank, and a load vector that will 
  // contain the corresponding loads. 
  int64_t* partition_vector = polymec_malloc(sizeof(int64_t) * orig_array_size);

  // Invert the send map to tell us the destination process for each index.
  int_int_unordered_map_t* proc_for_index = int_int_unordered_map_new();
  {
    int pos = 0, proc, *indices;
    while (int_ptr_unordered_map_next(sends, &pos, &proc, (void**)&indices))
    {
      for (int i = 0; i < num_points_to_rank[proc]; ++i)
        int_int_unordered_map_insert(proc_for_index, indices[i], proc);
    }
  }
  int_ptr_unordered_map_free(sends);

  // Now fill the partition vector.
  for (int i = 0; i < orig_array_size; ++i)
  {
    int* proc_ptr = int_int_unordered_map_get(proc_for_index, i);
    if (proc_ptr != NULL)
    {
      partition_vector[i] = (int64_t)(*proc_ptr);
    }
    else
      partition_vector[i] = (int64_t)rank;
  }

  // Clean up.
  int_int_unordered_map_free(proc_for_index);

  STOP_FUNCTION_TIMER();
  return partition_vector;
}

#endif

int64_t* repartition_points(point_t* local_points, 
                            size_t num_local_points,
                            MPI_Comm comm,
                            int* weights,
                            real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);

#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
  {
    int64_t* P = polymec_calloc(num_local_points, sizeof(int64_t));
    STOP_FUNCTION_TIMER();
    return P;
  }

  // Set up a Hilbert space filling curve that can map the given points to indices.
  hilbert_t* hilbert = hilbert_from_points(local_points, num_local_points);

  // Create an array of 4-tuples containing the 
  // (Hilbert index, rank, index, weight) of each point. 
  // Partitioning the points amounts to sorting this array and breaking it into 
  // parts whose work is equal. 
  index_t* part_array = polymec_malloc(sizeof(index_t) * 4 * num_local_points);
  for (int i = 0; i < num_local_points; ++i)
  {
    part_array[4*i] = hilbert_index(hilbert, &local_points[i]);
    part_array[4*i+1] = (index_t)rank;
    part_array[4*i+2] = (index_t)i;
    part_array[4*i+3] = 1;
  }
  if (weights != NULL)
  {
    for (int i = 0; i < num_local_points; ++i)
      part_array[4*i+3] = (index_t)weights[i];
  }

  // Now we have a distributed array, stored in segments on the processors 
  // in this communicator. Sort it so that process p holds the elements (in 
  // ascending order) that are greater than those of p-1 and less than those 
  // of p+1.
  parallel_sort(comm, part_array, num_local_points, 
                4*sizeof(index_t), hilbert_comp);
  release_ref(hilbert);

  // Try to balance the workload and bug out if we fail.
  size_t balanced_num_points = num_local_points;
  bool balanced = balance_loads(comm, imbalance_tol, &part_array, &balanced_num_points);
  if (!balanced)
  {
    polymec_free(part_array);
    return NULL;
  }

  // Now we create local partition/load vectors for each process using the elements
  // in the sorted list.
  int64_t* local_partition = partition_vector_from_balanced_array(comm, part_array, 
                                                                  balanced_num_points,
                                                                  num_local_points);

  // Clean up.
  polymec_free(part_array);

  STOP_FUNCTION_TIMER();
  return local_partition;
#else
  int64_t* P = polymec_calloc(num_local_points, sizeof(int64_t));
  return P;
#endif
}

redistribution_t* redistribution_from_partition(MPI_Comm comm, 
                                                int64_t* local_partition,
                                                size_t num_local_vertices)
{
  // Tally up the number of vertices we're going to send to every other process, 
  // including those that are staying on our own.
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  int num_vertices_to_send[nprocs];
  memset(num_vertices_to_send, 0, sizeof(int) * nprocs);
  for (int v = 0; v < (int)num_local_vertices; ++v)
    num_vertices_to_send[local_partition[v]]++;

  int num_vertices_to_receive[nprocs];
  MPI_Alltoall(num_vertices_to_send, 1, MPI_INT, 
               num_vertices_to_receive, 1, MPI_INT, comm);

  // Set up send and receive vertices.
  int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* receive_map = int_ptr_unordered_map_new();
  int num_local_vert = num_vertices_to_receive[rank];
  int recv_offset = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    if (p != rank)
    {
      if ((num_vertices_to_receive[p] > 0) && (p != rank))
      {
        if (p != rank)
        {
          int num_vert = num_vertices_to_receive[p];
          int_array_t* receive_vertices = int_array_new_with_size(num_vertices_to_receive[p]);
          for (int i = 0; i < num_vert; ++i, ++recv_offset)
            receive_vertices->data[i] = num_local_vert + recv_offset;
          int_ptr_unordered_map_insert(receive_map, p, receive_vertices);
        }
      }
      if (num_vertices_to_send[p] > 0)
      {
        int_array_t* send_vertices = int_array_new();
        for (int v = 0; v < num_local_vertices; ++v)
        {
          if (local_partition[v] == p)
            int_array_append(send_vertices, v);
        }
        if (p != rank)
          int_ptr_unordered_map_insert(send_map, p, send_vertices);
      }
    }
  }

  // Stick everything into a redistribution object.
  redistribution_t* redist = polymec_malloc(sizeof(redistribution_t));
  redist->send_procs = int_array_new_with_size(send_map->size);
  redist->send_indices = polymec_malloc(sizeof(int_array_t*) * send_map->size);
  int pos = 0, proc, i = 0;
  void* val;
  while (int_ptr_unordered_map_next(send_map, &pos, &proc, &val))
  {
    redist->send_procs->data[i] = proc;
    redist->send_indices[i] = (int_array_t*)val;
    ++i;
  }
  int_ptr_unordered_map_free(send_map);

  redist->receive_procs = int_array_new_with_size(receive_map->size);
  redist->receive_indices = polymec_malloc(sizeof(int_array_t*) * receive_map->size);
  pos = 0; 
  i = 0;
  while (int_ptr_unordered_map_next(receive_map, &pos, &proc, &val))
  {
    redist->receive_procs->data[i] = proc;
    redist->receive_indices[i] = (int_array_t*)val;
    ++i;
  }
  int_ptr_unordered_map_free(receive_map);
  return redist;
}

void redistribution_free(redistribution_t* redist)
{
  for (size_t p = 0; p < redist->send_procs->size; ++p)
    int_array_free(redist->send_indices[p]);
  polymec_free(redist->send_indices);
  int_array_free(redist->send_procs);
  for (size_t p = 0; p < redist->receive_procs->size; ++p)
    int_array_free(redist->receive_indices[p]);
  polymec_free(redist->receive_indices);
  int_array_free(redist->receive_procs);
  polymec_free(redist);
}

