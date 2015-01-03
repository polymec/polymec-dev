// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_point_cloud.h"
#include "core/hilbert.h"
#include "core/parallel_qsort.h"

#if POLYMEC_HAVE_MPI

// This helper constructs and returns a point cloud from the points with the 
// given indices in the given point cloud.
static point_cloud_t* create_subcloud(MPI_Comm comm, 
                                      point_cloud_t* cloud, 
                                      int* indices, int num_indices)
{
  // This is super easy--just pick out the points we want!
  point_cloud_t* subcloud = point_cloud_new(comm, num_indices);
  for (int i = 0; i < num_indices; ++i)
    subcloud->points[i] = cloud->points[indices[i]];
  return subcloud;
}

static void point_cloud_migrate(point_cloud_t** cloud, 
                                exchanger_t* migrator)
{
}

// This helper creates a Hilbert space-filling curve from a set of points.
static hilbert_t* hilbert_from_points(point_t* points, int num_points)
{
  bbox_t bbox = {.x1 = FLT_MAX, .x2 = -FLT_MAX, 
                 .y1 = FLT_MAX, .y2 = -FLT_MAX, 
                 .z1 = FLT_MAX, .z2 = -FLT_MAX};
  for (int i = 0; i < num_points; ++i)
    bbox_grow(&bbox, &points[i]);

  // (Handle lower dimensional point distributions gracefully.)
  if (fabs(bbox.x2 - bbox.x1) < FLT_MIN)
  {
    bbox.x1 -= 0.5;
    bbox.x2 += 0.5;
  }
  if (fabs(bbox.y2 - bbox.y1) < FLT_MIN)
  {
    bbox.y1 -= 0.5;
    bbox.y2 += 0.5;
  }
  if (fabs(bbox.z2 - bbox.z1) < FLT_MIN)
  {
    bbox.z1 -= 0.5;
    bbox.z2 += 0.5;
  }
  return hilbert_new(&bbox);
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

// This creates a local partition vector using the information in the sorted 
// distributed array. This is to be used only with repartition_point_cloud().
static int* create_partition_from_sorted_array(MPI_Comm comm, 
                                               index_t* array, 
                                               int local_array_size)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Find out who expects data from us.
  int num_points_for_rank[nprocs], my_num_points_for_rank[nprocs];
  memset(my_num_points_for_rank, 0, sizeof(int) * nprocs);
  for (int i = 0; i < local_array_size; ++i)
  {
    int rank = array[4*i+1];
    ++my_num_points_for_rank[rank];
  }
  MPI_Alltoall(my_num_points_for_rank, 1, MPI_INT, num_points_for_rank, 1, MPI_INT, comm);

  // Now find out which ranks got our points.
  // FIXME
  return NULL;
}

#endif

exchanger_t* distribute_point_cloud(point_cloud_t** cloud, 
                                    MPI_Comm comm,
                                    int64_t* global_partition)
{
#if POLYMEC_HAVE_MPI
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Make sure we're all here.
  MPI_Barrier(comm);

  point_cloud_t* global_cloud = *cloud;
  int num_global_points = (rank == 0) ? global_cloud->num_points : 0;
  point_cloud_t* local_cloud = NULL;
  uint64_t vtx_dist[nprocs+1];
  if (rank == 0)
  {
    // Take stock of how many points we'll have per process.
    int num_points[nprocs];
    memset(num_points, 0, sizeof(int) * nprocs);
    for (int i = 0; i < global_cloud->num_points; ++i)
      num_points[global_partition[i]]++;

    // Construct the distribution of vertices for the partitioning.
    vtx_dist[0] = 0;
    for (int p = 0; p < nprocs; ++p)
      vtx_dist[p+1] = vtx_dist[p] + num_points[p];

    // Carve out the portion of the cloud that will stick around on process 0.
    {
      int indices[num_points[0]], k = 0;
      for (int i = 0; i < global_cloud->num_points; ++i)
      {
        if (global_partition[i] == rank)
          indices[k++] = i;
      }
      local_cloud = create_subcloud(comm, global_cloud, indices, num_points[0]);
    }

    // Now do the other processes.
    serializer_t* ser = point_cloud_serializer();
    byte_array_t* bytes = byte_array_new();
    for (int p = 1; p < nprocs; ++p)
    {
      // Share the vertex distribution.
      MPI_Send(vtx_dist, nprocs+1, MPI_UINT64_T, p, p, comm);

      // Create the pth subcloud.
      int indices[num_points[p]], k = 0;
      for (int i = 0; i < global_cloud->num_points; ++i)
      {
        if (global_partition[i] == p)
          indices[k++] = i;
      }
      point_cloud_t* p_cloud = create_subcloud(comm, global_cloud, indices, num_points[p]);

      // Serialize it and send its size (and it) to process p.
      size_t offset = 0;
      serializer_write(ser, p_cloud, bytes, &offset);
      MPI_Send(&bytes->size, 1, MPI_INT, p, p, comm);
      MPI_Send(bytes->data, bytes->size, MPI_BYTE, p, p, comm);

      // Clean up.
      byte_array_clear(bytes);
      point_cloud_free(p_cloud);
    }
    ser = NULL;
    byte_array_free(bytes);
  }
  else
  {
    // Receive the vertex distribution of the incoming cloud.
    MPI_Status status;
    MPI_Recv(vtx_dist, nprocs+1, MPI_UINT64_T, 0, rank, comm, &status);

    // Receive the size of the incoming cloud.
    int cloud_size;
    MPI_Recv(&cloud_size, 1, MPI_INT, 0, rank, comm, &status);

    // Now receive the cloud.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, cloud_size);

    MPI_Recv(bytes->data, cloud_size, MPI_BYTE, 0, rank, comm, &status);
    serializer_t* ser = point_cloud_serializer();
    size_t offset = 0;
    local_cloud = serializer_read(ser, bytes, &offset);
    
    byte_array_free(bytes);
    ser = NULL;
  }

  *cloud = local_cloud;

  // Clean up.
  if (global_cloud != NULL)
    point_cloud_free(global_cloud);

  return create_distributor(comm, global_partition, num_global_points);
#else
  return exchanger_new(comm);
#endif
}

int64_t* partition_vector_from_point_cloud(point_cloud_t* global_cloud, 
                                           MPI_Comm comm, 
                                           int* weights, 
                                           real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_cloud != NULL));

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
  {
    // Dumb, but correct.
    int64_t* global_partition = polymec_malloc(sizeof(int64_t) * global_cloud->num_points);
    memset(global_partition, 0, sizeof(int64_t) * global_cloud->num_points);
    return global_partition;
  }

#ifndef NDEBUG
  // Make sure there are enough points to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT(global_cloud->num_points > nprocs);
  }
#endif

  int64_t* global_partition = NULL;
  int num_global_points = 0;
  if (rank == 0)
  {
    num_global_points = global_cloud->num_points;

    // Set up a Hilbert space filling curve that can map the given points to indices.
    hilbert_t* hilbert = hilbert_from_points(global_cloud->points, num_global_points);

    // Create an array of 2-tuples containing the (Hilbert index, weight) of each point. 
    // Partitioning the points amounts to sorting this array and breaking it into parts whose 
    // work is equal. Also sum up the work on the points.
    index_t* part_array = polymec_malloc(sizeof(index_t) * 3 * num_global_points);
    uint64_t total_work = 0;
    if (weights != NULL)
    {
      for (int i = 0; i < num_global_points; ++i)
      {
        part_array[3*i]   = i;
        part_array[3*i+1] = hilbert_index(hilbert, &global_cloud->points[i]);
        part_array[3*i+2] = (index_t)weights[i];
        total_work += weights[i];
      }
    }
    else
    {
      for (int i = 0; i < num_global_points; ++i)
      {
        part_array[3*i]   = i;
        part_array[3*i+1] = hilbert_index(hilbert, &global_cloud->points[i]);
        part_array[3*i+2] = 1;
      }
      total_work = num_global_points;
    }

    // Sort the array.
    qsort(part_array, (size_t)num_global_points, 3*sizeof(index_t), hilbert_comp);

    // Now we need to break it into parts of equal work.
    real_t work_per_proc = 1.0 * total_work / nprocs;
    int part_offsets[nprocs+1];
    real_t part_work[nprocs+1];
    part_offsets[0] = 0;
    part_work[0] = 0.0;
    for (int p = 0; p < nprocs; ++p)
    {
      int i = part_offsets[p];
      real_t work = 0.0, last_weight = 0.0, cum_work = part_work[p];
      while ((cum_work < ((p+1) * work_per_proc)) && 
             (i < num_global_points))
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
    global_partition = polymec_malloc(sizeof(int64_t) * num_global_points);
    int k = 0;
    for (int p = 0; p < nprocs; ++p)
    {
      for (int i = part_offsets[p]; i < part_offsets[p+1]; ++i, ++k)
      {
        index_t j = part_array[3*k];
        global_partition[j] = p;
      }
    }
    
    // Clean up.
    polymec_free(part_array);
  }

  return global_partition;

#else
  // This is dumb, but we were asked for it.
  int64_t* global_partition = polymec_malloc(sizeof(int64_t) * global_cloud->num_points);
  memset(global_partition, 0, sizeof(int64_t) * global_cloud->num_points);
  return global_partition;
#endif
}

exchanger_t* partition_point_cloud(point_cloud_t** cloud, MPI_Comm comm, int* weights, real_t imbalance_tol)
{
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  point_cloud_t* cl = *cloud;

#if POLYMEC_HAVE_MPI
  ASSERT((*cloud == NULL) || ((*cloud)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  // Now do the space-filling curve thing.
  int64_t* global_partition = partition_vector_from_point_cloud(cl, comm, weights, imbalance_tol);

  // Distribute the point cloud.
  exchanger_t* distributor = distribute_point_cloud(cloud, comm, global_partition);

  // Clean up.
  polymec_free(global_partition);

  // Return the migrator.
  return (distributor == NULL) ? exchanger_new(comm) : distributor;
#else
  return exchanger_new(cl->comm);
#endif
}

#if 0
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

  // Set up a Hilbert space filling curve that can map the given points to indices.
  hilbert_t* hilbert = hilbert_from_points(cl->points, cl->num_points);

  // Create an array of 4-tuples containing the 
  // (Hilbert index, rank, index, weight) of each point. 
  // Partitioning the points amounts to sorting this array and breaking it into parts whose 
  // work is equal. Also sum up the work on the points.
  index_t* part_array = polymec_malloc(sizeof(index_t) * 4 * cl->num_points);
  uint64_t total_work = 0;
  if (weights != NULL)
  {
    for (int i = 0; i < cl->num_points; ++i)
    {
      part_array[4*i] = hilbert_index(hilbert, &cl->points[i]);
      part_array[4*i+1] = (index_t)rank;
      part_array[4*i+2] = (index_t)i;
      part_array[4*i+3] = (index_t)weights[i];
      total_work += weights[i];
    }
  }
  else
  {
    for (int i = 0; i < cl->num_points; ++i)
    {
      part_array[4*i] = hilbert_index(hilbert, &cl->points[i]);
      part_array[4*i+1] = (index_t)rank;
      part_array[4*i+2] = (index_t)i;
      part_array[4*i+1] = 1;
    }
    total_work = cl->num_points;
  }

  // Now we have a distributed array, stored in segments on the processors in this communicator.
  // Sort the thing all-parallel-like using regular sampling.
  parallel_qsort(cl->comm, part_array, cl->num_points, 
                 4*sizeof(index_t), hilbert_comp, NULL);

  // Now we create a local partition vector for each process using the elements
  // in the sorted list.
  int* local_partition = create_partition_from_sorted_array(cl->comm, part_array, cl->num_points);

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
#endif

