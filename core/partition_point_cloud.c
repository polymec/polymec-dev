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

#include "core/partition_point_cloud.h"
#include "core/hilbert.h"

// This helper constructs and returns a point cloud from the given points in 
// the given point cloud.
static point_cloud_t* create_subcloud(point_cloud_t* cloud, int* indices, int num_indices)
{
  return NULL;
}

static void point_cloud_distribute(point_cloud_t** cloud, 
                                   MPI_Comm comm,
                                   int64_t* global_partition)
{
}

static void point_cloud_migrate(point_cloud_t** cloud, 
                                exchanger_t* migrator)
{
}

#if POLYMEC_HAVE_MPI
// This helper is a comparison function used to sort (Hilbert index, weight) tuples.
// Only the Hilbert index factors into the ordering.
static int hilbert_comp(const void* l, const void* r)
{
  const index_t* li = l;
  const index_t* ri = r;
  return (li[0] < ri[0]) ? -1
                         : (li[0] > ri[0]) ? 1
                                           : 0;
}
#endif

exchanger_t* partition_point_cloud(point_cloud_t** cloud, MPI_Comm comm, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  point_cloud_t* cl = *cloud;

#if POLYMEC_HAVE_MPI

  int nprocs, rank;
  MPI_Comm_size(cl->comm, &nprocs);
  MPI_Comm_rank(cl->comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(cl->comm);

  // Clouds on rank != 0 must be NULL.
  ASSERT((rank == 0) || (*cloud == NULL));

  int64_t* global_partition = NULL;
  if (rank == 0)
  {
    // Set up a Hilbert space filling curve that can map the given points to indices.
    // Also, sum up all the work on the points.
    bbox_t bbox = {.x1 = FLT_MAX, .x2 = -FLT_MAX, .y1 = FLT_MAX, .y2 = -FLT_MAX, 
      .z1 = FLT_MAX, .z2 = -FLT_MAX};
    for (int i = 0; i < cl->num_points; ++i)
      bbox_grow(&bbox, &cl->points[i]);
    hilbert_t* hilbert = hilbert_new(&bbox);

    // Create an array of 2-tuples containing the (Hilbert index, weight) of each point. 
    // Partitioning the points amounts to sorting this array and breaking it into parts whose 
    // work is equal. Also sum up the work on the points.
    index_t* part_array = polymec_malloc(sizeof(index_t) * 2 * cl->num_points);
    uint64_t total_work = 0;
    if (weights != NULL)
    {
      for (int i = 0; i < cl->num_points; ++i)
      {
        part_array[2*i] = hilbert_index(hilbert, &cl->points[i]);
        part_array[2*i+1] = (index_t)weights[i];
        total_work += weights[i];
      }
    }
    else
    {
      for (int i = 0; i < cl->num_points; ++i)
      {
        part_array[2*i] = hilbert_index(hilbert, &cl->points[i]);
        part_array[2*i+1] = 1;
      }
      total_work = cl->num_points;
    }

    // Sort the array.
    qsort(part_array, (size_t)cl->num_points, 2*sizeof(index_t), hilbert_comp);

    // Now we need to break it into parts of equal work.
    real_t work_per_proc = 1.0 * total_work / nprocs;
    int part_offsets[nprocs+1];
    part_offsets[0] = 0;
    for (int p = 0; p < nprocs; ++p)
    {
      int i = part_offsets[p];
      real_t work = 0.0;
      do
      {
        work += 1.0 * part_array[2*i+1];
        ++i;
      }
      while (work < work_per_proc);

      // If we've obviously overloaded this process, back up one step.
      if (((work - work_per_proc)/work_per_proc > imbalance_tol) && 
          ((work_per_proc - (work - weights[i-1]) <= imbalance_tol)))
        --i;
      part_offsets[p+1] = i;
    }
    
    // Now we create the global partition vector and fill it.
    global_partition = polymec_malloc(sizeof(int64_t) * cl->num_points);
    int k = 0;
    for (int p = 0; p < nprocs; ++p)
    {
      for (int i = part_offsets[p]; i < part_offsets[p+1]; ++i, ++k)
        global_partition[k] = p;
    }
  }

  // Distribute the point cloud.
  point_cloud_distribute(cloud, comm, global_partition);

  // Set up an exchanger to distribute field data.
  int num_points = (cl != NULL) ? cl->num_points : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_points);

  // Clean up.
  polymec_free(global_partition);

  // Return the migrator.
  return (distributor == NULL) ? exchanger_new(cl->comm) : distributor;
#else
  return exchanger_new(cl->comm);
#endif
}

exchanger_t* repartition_point_cloud(point_cloud_t** cloud, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  point_cloud_t* cl = *cloud;

#if POLYMEC_HAVE_MPI
  int nprocs, rank;
  MPI_Comm_size(cl->comm, &nprocs);
  MPI_Comm_rank(cl->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(cl->comm);

#if 0
  // Map the graph to the different domains, producing a local partition vector
  // (with values included for ghost cells).
  int* local_partition = repartition_graph(local_graph, cl->num_ghost_points, 
                                                  weights, imbalance_tol, cloud_ex);

  // Set up an exchanger to migrate field data.
  int num_vertices = adj_graph_num_vertices(local_graph);
  exchanger_t* migrator = create_migrator(cl->comm, local_partition, num_vertices);

  // Migrate the point cloud.
  point_cloud_migrate(cloud, local_graph, migrator);

  // Clean up.
  adj_graph_free(local_graph);
  polymec_free(local_partition);

  // Return the migrator.
  return (migrator == NULL) ? exchanger_new(cl->comm) : migrator;
#endif
  return NULL;
#else
  return exchanger_new(cl->comm);
#endif
}

