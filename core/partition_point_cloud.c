// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_point_cloud.h"
#include "core/hilbert.h"
#include "core/parallel_sort.h"
#include "core/unordered_set.h"
#include "core/timer.h"

#if POLYMEC_HAVE_MPI

// This helper constructs and returns a point cloud from the points with the 
// given indices in the given point cloud.
static point_cloud_t* create_subcloud(MPI_Comm comm, 
                                      point_cloud_t* cloud, 
                                      int* indices, int num_indices)
{
  START_FUNCTION_TIMER();
  // This is super easy--just pick out the points we want!
  point_cloud_t* subcloud = point_cloud_new(comm, num_indices);
  for (int i = 0; i < num_indices; ++i)
    subcloud->points[i] = cloud->points[indices[i]];

  // Now copy properties using their serializers.
  int pos = 0;
  char* prop_name;
  void* prop_data;
  serializer_t* prop_ser;
  while (point_cloud_next_property(cloud, &pos, &prop_name, &prop_data, &prop_ser))
  {
    if (prop_ser != NULL)
    {
      void* prop_data_copy = serializer_clone_object(prop_ser, prop_data);
      point_cloud_set_property(subcloud, prop_name, prop_data_copy, prop_ser);
    }
    else
      log_debug("create_subcloud: property '%s' has no serializer.", prop_name);
  }

  STOP_FUNCTION_TIMER();
  return subcloud;
}

// Fuse a set of subclouds into a single point cloud. Ghost points are not 
// permitted.
static point_cloud_t* fuse_clouds(point_cloud_t** subclouds, int num_subclouds)
{
  int num_points = 0;
  for (int i = 0; i < num_subclouds; ++i)
  {
    ASSERT(subclouds[i]->num_ghosts == 0);
    num_points += subclouds[i]->num_points;
  }

  point_cloud_t* fused_cloud = point_cloud_new(subclouds[0]->comm, num_points);
  int k = 0;
  for (int i = 0; i < num_subclouds; ++i)
  {
    for (int j = 0; j < subclouds[i]->num_points; ++j, ++k)
      fused_cloud->points[k] = subclouds[i]->points[j];
  }

  return fused_cloud;
}

// Migrate point cloud data using the given exchanger.
static void point_cloud_migrate(point_cloud_t** cloud, 
                                exchanger_t* migrator)
{
  START_FUNCTION_TIMER();
  point_cloud_t* c = *cloud;

  // Post receives for buffer sizes.
  int num_receives = exchanger_num_receives(migrator);
  int num_sends = exchanger_num_sends(migrator);
  int receive_buffer_sizes[num_receives], receive_procs[num_receives];
  int pos = 0, proc, num_indices, *indices, i_req = 0;
  MPI_Request requests[num_receives + num_sends];
  while (exchanger_next_receive(migrator, &pos, &proc, &indices, &num_indices))
  {
    receive_procs[i_req] = proc;
    MPI_Irecv(&receive_buffer_sizes[i_req], 1, MPI_INT, proc, 0, c->comm, &requests[i_req]);
    ++i_req;
  }

  // Build point clouds to send to other processes.
  int_unordered_set_t* sent_points = int_unordered_set_new();
  serializer_t* ser = point_cloud_serializer();
  byte_array_t* send_buffers[num_sends];
  int send_procs[num_sends];
  pos = 0;
  while (exchanger_next_send(migrator, &pos, &proc, &indices, &num_indices))
  {
    send_procs[i_req-num_receives] = proc;
    byte_array_t* bytes = byte_array_new();

    // Add the indices of the cells we are sending.
    for (int i = 0; i < num_indices; ++i)
      int_unordered_set_insert(sent_points, indices[i]);

    // Create the subcloud to send. 
    point_cloud_t* subcloud = create_subcloud(c->comm, c, indices, num_indices);

    // Serialize and send the buffer size.
    size_t offset = 0;
    serializer_write(ser, subcloud, bytes, &offset);
    MPI_Isend(&bytes->size, 1, MPI_INT, proc, 0, c->comm, &requests[i_req]);

    // Clean up.
    point_cloud_free(subcloud);
    send_buffers[i_req - num_receives] = bytes;
    ++i_req;
  }
  ASSERT(i_req == num_sends + num_receives);

  // Wait for the buffer sizes to be transmitted.
  MPI_Status statuses[num_receives + num_sends];
  MPI_Waitall(num_receives + num_sends, requests, statuses);

  // Post receives for the actual messages.
  byte_array_t* receive_buffers[num_receives];
  for (int i = 0; i < num_receives; ++i)
  {
    receive_buffers[i] = byte_array_new();
    byte_array_resize(receive_buffers[i], receive_buffer_sizes[i]);
    MPI_Irecv(receive_buffers[i]->data, receive_buffer_sizes[i], MPI_BYTE, receive_procs[i], 0, c->comm, &requests[i]);
  }

  // Send the actual clouds and wait for receipt.
  for (int i = 0; i < num_sends; ++i)
    MPI_Isend(send_buffers[i]->data, send_buffers[i]->size, MPI_BYTE, send_procs[i], 0, c->comm, &requests[num_receives + i]);
  MPI_Waitall(num_receives + num_sends, requests, statuses);

  // Unpack the clouds.
  point_cloud_t* subclouds[1+num_receives];
  for (int i = 0; i < num_receives; ++i)
  {
    size_t offset = 0;
    subclouds[i+1] = serializer_read(ser, receive_buffers[i], &offset);
  }

  // Clean up all the stuff from the exchange.
  ser = NULL;
  for (int i = 0; i < num_receives; ++i)
    byte_array_free(receive_buffers[i]);
  for (int i = 0; i < num_sends; ++i)
    byte_array_free(send_buffers[i]);

  // Construct a local subcloud and store it in subclouds[0]. This subcloud
  // consists of all points not sent to other processes.
  {
    int num_local_points = c->num_points - sent_points->size;
    int local_points[num_local_points], j = 0;
    for (int i = 0; i < c->num_points; ++i)
    {
      if (!int_unordered_set_contains(sent_points, i))
        local_points[j++] = i;
    }
    subclouds[0] = create_subcloud(c->comm, c, local_points, num_local_points);
  }

  // Fuse all the subclouds into a single point cloud.
  int_unordered_set_free(sent_points);
  point_cloud_free(c);
  *cloud = fuse_clouds(subclouds, 1+num_receives);
  STOP_FUNCTION_TIMER();
  POLYMEC_NOT_IMPLEMENTED
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

#endif

exchanger_t* distribute_point_cloud(point_cloud_t** cloud, 
                                    MPI_Comm comm,
                                    int64_t* global_partition)
{
#if POLYMEC_HAVE_MPI
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  // Make sure we're all here.
  MPI_Barrier(comm);

  START_FUNCTION_TIMER();

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

  exchanger_t* ex = create_distributor(comm, global_partition, num_global_points);
  STOP_FUNCTION_TIMER();
  return ex;
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
  START_FUNCTION_TIMER();
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
    STOP_FUNCTION_TIMER();
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

    // Create an array of 3-tuples containing the 
    // (point index, Hilbert index, weight) of each point. 
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

  STOP_FUNCTION_TIMER();
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
  START_FUNCTION_TIMER();
  ASSERT((*cloud == NULL) || ((*cloud)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
  {
    STOP_FUNCTION_TIMER();
    return exchanger_new(comm);
  }

  // Now do the space-filling curve thing.
  int64_t* global_partition = partition_vector_from_point_cloud(cl, comm, weights, imbalance_tol);

  // Distribute the point cloud.
  exchanger_t* distributor = distribute_point_cloud(cloud, comm, global_partition);

  // Clean up.
  polymec_free(global_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return (distributor == NULL) ? exchanger_new(comm) : distributor;
#else
  return exchanger_new(cl->comm);
#endif
}

//------------------------------------------------------------------------
//                     Dynamic repartition code below
//------------------------------------------------------------------------

#if POLYMEC_HAVE_MPI 

// This creates local partition and load vectors using the information in the sorted 
// distributed array. This is to be used only with repartition_point_cloud().
static void create_partition_and_load_from_sorted_array(MPI_Comm comm, 
                                                        index_t* array, 
                                                        int local_array_size,
                                                        int64_t** partition_vector,
                                                        index_t** load_vector)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Find out who is sending us data by reading off the process ranks in our local 
  // portion of the sorted array.

  // Store the number of points we are receiving from rank p in 
  // num_points_from_rank[p].
  int num_points_from_rank[nprocs];
  {
    int my_num_points_from_rank[nprocs];
    memset(my_num_points_from_rank, 0, sizeof(int) * nprocs);
    for (int i = 0; i < local_array_size; ++i)
    {
      int rank = array[4*i+1];
      ++my_num_points_from_rank[rank];
    }
    MPI_Alltoall(my_num_points_from_rank, 1, MPI_INT, num_points_from_rank, 1, MPI_INT, comm);
  }

  // Store the number of points we are sending to rank p in 
  // num_points_to_rank[p].
  int num_points_to_rank[nprocs];
  MPI_Alltoall(num_points_from_rank, 1, MPI_INT, num_points_to_rank, 1, MPI_INT, comm);
  // FIXME: Does this work???

  // Post receives/sends for all nonzero sets of points.
  for (int p = 0; p < nprocs; ++p)
  {
    if (num_points_from_rank[p] > 0)
    {
      // FIXME
    }

    if (num_points_to_rank[p] > 0)
    {
      // FIXME
    }
  }

  // Wait for them to send.

  // FIXME
  POLYMEC_NOT_IMPLEMENTED
  STOP_FUNCTION_TIMER();
  return NULL;
}

// This rebalances the workload on the given set of points, changing the distribution 
// of the points on the processors but preserving their ordering. The points are 
// assumed to be 4-wide sets of indices as ordered within repartition_point_cloud
// below. This procedure is serial, so it must be used sparingly.
static void balance_load(MPI_Comm comm, 
                         int64_t* local_partition, 
                         index_t* local_load,
                         int num_local_points,
                         real_t imbalance_tol)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Sum the total load on all processes, and find the ideal load per process.
  real_t my_load = 0.0;
  for (int i = 0; i < num_local_points; ++i)
    my_load += 1.0 * local_load[i];
  real_t total_load;
  MPI_Allreduce(&my_load, &total_load, 1, MPI_REAL_T, MPI_SUM, comm);
  real_t ideal_proc_load = total_load / nprocs;

  // Now loop through the processes and carve off the ideal workload, .
  for (int p = 0; p < nprocs; ++p)
  {
    if (rank == p)
    {

    }
  }
}

#endif

exchanger_t* repartition_point_cloud(point_cloud_t** cloud, 
                                     int* weights, 
                                     real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  point_cloud_t* cl = *cloud;

#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(cl->comm, &nprocs);
  MPI_Comm_rank(cl->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
  {
    STOP_FUNCTION_TIMER();
    return exchanger_new(cl->comm);
  }

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
      part_array[4*i+3] = 1;
    }
    total_work = cl->num_points;
  }

  // Now we have a distributed array, stored in segments on the processors 
  // in this communicator. Sort it so that process p holds the elements (in 
  // ascending order) that are greater than those of p-1 and less than those 
  // of p+1.
  parallel_sort(cl->comm, part_array, cl->num_points, 
                4*sizeof(index_t), hilbert_comp);

  // Now we create local partition/load vectors for each process using the elements
  // in the sorted list.
  int64_t* local_partition;
  index_t* local_load;
  create_partition_and_load_from_sorted_array(cl->comm, part_array, cl->num_points,
                                             &local_partition, &local_load);

  if (weights != NULL)
  {
    // We make adjustments to the partition vector to accomodate variable loads.
    balance_load(cl->comm, local_partition, local_load, 
                 cl->num_points, imbalance_tol);
  }

  // Set up an exchanger to migrate field data.
  exchanger_t* migrator = create_migrator(cl->comm, local_partition, cl->num_points);

  // Migrate the point cloud.
  point_cloud_migrate(cloud, migrator);

  // Clean up.
  polymec_free(local_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return (migrator == NULL) ? exchanger_new(cl->comm) : migrator;
#else
  return exchanger_new(cl->comm);
#endif
}

