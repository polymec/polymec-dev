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

migrator_t* distribute_point_cloud(point_cloud_t** cloud, 
                                   MPI_Comm comm,
                                   int64_t* global_partition)
{
#if POLYMEC_HAVE_MPI
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return migrator_new(comm);

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

  migrator_t* m = migrator_from_global_partition(comm, global_partition, num_global_points);
  STOP_FUNCTION_TIMER();
  return m;
#else
  return migrator_new(comm);
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

migrator_t* partition_point_cloud(point_cloud_t** cloud, MPI_Comm comm, int* weights, real_t imbalance_tol)
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
    return migrator_new(comm);
  }

  // Now do the space-filling curve thing.
  int64_t* global_partition = partition_vector_from_point_cloud(cl, comm, weights, imbalance_tol);

  // Distribute the point cloud.
  migrator_t* migrator = distribute_point_cloud(cloud, comm, global_partition);

  // Clean up.
  polymec_free(global_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return (migrator == NULL) ? migrator_new(comm) : migrator;
#else
  return migrator_new(cl->comm);
#endif
}

//------------------------------------------------------------------------
//                     Dynamic repartition code below
//------------------------------------------------------------------------

#if POLYMEC_HAVE_MPI 

// This rebalances the workload on the given set of points, changing the distribution 
// of the points on the processors but preserving their ordering. The points are 
// assumed to be 4-wide sets of indices as ordered within repartition_point_cloud
// below. This procedure is serial, so the function must be used sparingly. 
// It returns true if it can balance the loads to within the imbalance toleranace,
// false if not.
static bool balance_loads(MPI_Comm comm, 
                          real_t imbalance_tol,
                          index_t** local_array,
                          int* local_array_size)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  index_t* array = *local_array;
  int array_size = *local_array_size;

  // Sum the total load on all processes, and find the ideal load per process.
  real_t my_load = 0.0;
  for (int i = 0; i < array_size; ++i)
    my_load += 1.0 * array[4*i+3];
  log_debug("balance_loads: Current process workload is %g.", my_load);

  real_t total_load;
  MPI_Allreduce(&my_load, &total_load, 1, MPI_REAL_T, MPI_SUM, comm);
  int ideal_proc_load = (int)ceil(total_load / nprocs);
  log_debug("balance_loads: Ideal process workload is %d.", ideal_proc_load);

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
          int num_points_sent = 0, work_given = 0, delta;
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
      int load = 0;
      int num_points_kept = 0, delta;
      while ((load < ideal_proc_load) && (num_points_kept < array_size))
      {
        delta = (int)array[4*num_points_kept+3];
        load += delta;
        ++num_points_kept;
      }

      // Step back a point if it makes sense.
      if ((load > ideal_proc_load) && 
          (abs(load - delta - ideal_proc_load) < (load - ideal_proc_load)))
      {
        load -= delta;
        --num_points_kept;
      }

      // If there's a process "to our right," dump any excess work there.
      int excess_load = load - ideal_proc_load;
      if ((excess_load >= 0) || (num_points_kept < array_size))
      {
        if (rank == (nprocs - 1)) // last process!
        {
          // Can we fit the extra work on the process and still fall under 
          // our load imbalance tolerance?
          real_t imbalance = (1.0 * load - ideal_proc_load) / ideal_proc_load;
          if (imbalance > imbalance_tol)
          {
            log_debug("balance_loads: Load on last process (%d) is unbalanced (%g%%).", load, imbalance * 100.0);
            balanced = 0; 
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
          int num_points_sent = array_size - num_points_kept;
          ASSERT(num_points_sent >= 0);

          // Tell the sucker how many points we're giving him.
          log_debug("balance_loads: Sending %d points to next process.", num_points_sent);
          MPI_Send(&num_points_sent, 1, MPI_INT, rank+1, 0, comm);

          if (num_points_sent > 0)
          {
            // Pass the data along.
            MPI_Send(&array[4*num_points_kept], 4*num_points_sent, MPI_INDEX_T, 
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
      else if ((excess_load < 0) && (rank < (nprocs - 1))) // we need more work!
      {
        // We send the negation of how much work we would like to balance our 
        // load.
        int needed_work = (int)(ideal_proc_load - load);
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
        int last_load = (int)buffer[4*(num_points_received-1)+3];
        ASSERT(extra_work - last_load + load <= ideal_proc_load); // safety guarantee!
        if (extra_work + load > ideal_proc_load)
        {
          real_t imbalance = (1.0 * (load + excess_load) - ideal_proc_load) / ideal_proc_load;
          if (imbalance > imbalance_tol)
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
    }
  }

  // Let's all agree on whether the load was balanced.
  MPI_Allreduce(MPI_IN_PLACE, &balanced, 1, MPI_INT, MPI_MIN, comm);
  if (balanced == 1)
  {
    log_debug("balance_loads: Targeting %d local points.", array_size);
    log_debug("balance_loads: Successfully balanced workloads to within %d%%.", 
              (int)round(imbalance_tol * 100.0));
  }
  else
  {
    log_debug("balance_loads: Could not balance workloads to within %d%%.",
              (int)round(imbalance_tol * 100.0));
  }

  STOP_FUNCTION_TIMER();
  return (balanced == 1);
}

// This creates local partition and load vectors using the information in the sorted 
// distributed array. This is to be used only by repartition_point_cloud().
static int64_t* partition_vector_from_balanced_array(MPI_Comm comm, 
                                                     index_t* balanced_array, 
                                                     int balanced_array_size,
                                                     int orig_array_size)
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
          int rank = (int)balanced_array[4*i+1];
          if (rank == p)
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

// Migrate point cloud data using the given migrator.
static void point_cloud_migrate(point_cloud_t** cloud, 
                                migrator_t* migrator)
{
  START_FUNCTION_TIMER();
  point_cloud_t* c = *cloud;

  // FIXME: For now, we have to pull out the underlying exchanger for
  // FIXME: the migrator.
  exchanger_t* ex = migrator_exchanger(migrator);

  // Post receives for buffer sizes.
  int num_receives = exchanger_num_receives(ex);
  int num_sends = exchanger_num_sends(ex);
  int receive_buffer_sizes[num_receives], receive_procs[num_receives];
  int pos = 0, proc, num_indices, *indices, i_req = 0;
  MPI_Request requests[num_receives + num_sends];
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
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
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
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
}

#endif

migrator_t* repartition_point_cloud(point_cloud_t** cloud, 
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
    return migrator_new(cl->comm);
  }

  // Set up a Hilbert space filling curve that can map the given points to indices.
  hilbert_t* hilbert = hilbert_from_points(cl->points, cl->num_points);

  // Create an array of 4-tuples containing the 
  // (Hilbert index, rank, index, weight) of each point. 
  // Partitioning the points amounts to sorting this array and breaking it into 
  // parts whose work is equal. 
  index_t* part_array = polymec_malloc(sizeof(index_t) * 4 * cl->num_points);
  for (int i = 0; i < cl->num_points; ++i)
  {
    part_array[4*i] = hilbert_index(hilbert, &cl->points[i]);
    part_array[4*i+1] = (index_t)rank;
    part_array[4*i+2] = (index_t)i;
    part_array[4*i+3] = 1;
  }
  if (weights != NULL)
  {
    for (int i = 0; i < cl->num_points; ++i)
      part_array[4*i+3] = (index_t)weights[i];
  }

  // Now we have a distributed array, stored in segments on the processors 
  // in this communicator. Sort it so that process p holds the elements (in 
  // ascending order) that are greater than those of p-1 and less than those 
  // of p+1.
  parallel_sort(cl->comm, part_array, cl->num_points, 
                4*sizeof(index_t), hilbert_comp);

  // Try to balance the workload and bug out if we fail.
  int balanced_num_points = cl->num_points;
  bool balanced = balance_loads(cl->comm, imbalance_tol, &part_array, 
                                &balanced_num_points);
  if (!balanced)
  {
    polymec_free(part_array);
    return NULL;
  }

  // Now we create local partition/load vectors for each process using the elements
  // in the sorted list.
  int64_t* local_partition = partition_vector_from_balanced_array(cl->comm, part_array, 
                                                                  balanced_num_points,
                                                                  cl->num_points);

  // Set up an migrator to migrate field data.
  migrator_t* migrator = migrator_from_local_partition(cl->comm, local_partition, cl->num_points);

  // Migrate the point cloud.
  point_cloud_migrate(cloud, migrator);

  // Clean up.
  polymec_free(local_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return (migrator == NULL) ? migrator_new(cl->comm) : migrator;
#else
  return migrator_new(cl->comm);
#endif
}

