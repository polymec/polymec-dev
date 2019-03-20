// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partitioning.h"
#include "core/unordered_set.h"
#include "core/timer.h"
#include "geometry/partition_point_cloud.h"

#if POLYMEC_HAVE_MPI

// This helper constructs and returns a point cloud from the points with the
// given indices in the given point cloud.
static point_cloud_t* create_subcloud(MPI_Comm comm,
                                      point_cloud_t* cloud,
                                      int* indices, size_t num_indices)
{
  START_FUNCTION_TIMER();
  // This is super easy--just pick out the points we want!
  point_cloud_t* subcloud = point_cloud_new(comm, num_indices);
  for (int i = 0; i < num_indices; ++i)
    subcloud->points[i] = cloud->points[indices[i]];

  STOP_FUNCTION_TIMER();
  return subcloud;
}

#endif

void distribute_point_cloud(point_cloud_t** cloud,
                            MPI_Comm comm,
                            int64_t* global_partition,
                            point_cloud_field_t** fields,
                            size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning (or succeeds, if you like!).
  if (nprocs == 1)
    return;

  // Make sure we're all here.
  MPI_Barrier(comm);

  START_FUNCTION_TIMER();

  point_cloud_t* global_cloud = *cloud;
  point_cloud_t* local_cloud = NULL;
  point_cloud_field_t* local_fields[num_fields];
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

      // Handle fields.
      for (size_t i = 0; i < num_fields; ++i)
      {
        point_cloud_field_t* global_field = fields[i];
        point_cloud_field_t* local_field = point_cloud_field_new(local_cloud, global_field->num_components);
        int num_comps = local_field->num_components;
        for (size_t j = 0; j < local_field->num_local_values; ++j)
          for (size_t c = 0; c < num_comps; ++c)
            local_field->data[num_comps*j+c] = global_field->data[num_comps*indices[j]+c];
        local_fields[i] = local_field;
      }
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
      MPI_Send(bytes->data, (int)bytes->size, MPI_BYTE, p, p, comm);

      // Handle fields.
      for (size_t i = 0; i < num_fields; ++i)
      {
        point_cloud_field_t* global_field = fields[i];
        int num_comps = global_field->num_components;
        MPI_Send(&num_comps, 1, MPI_SIZE_T, p, p, comm);
        real_t local_vals[num_comps*num_points[p]];
        for (size_t j = 0; j < num_points[p]; ++j)
          for (size_t c = 0; c < num_comps; ++c)
            local_vals[num_comps*j+c] = global_field->data[num_comps*indices[j]+c];
        MPI_Send(local_vals, (int)(num_comps * num_points[p]), MPI_REAL_T, p, p, comm);
      }

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

    // Handle fields.
    for (size_t i = 0; i < num_fields; ++i)
    {
      int num_comps;
      MPI_Recv(&num_comps, 1, MPI_SIZE_T, 0, rank, comm, &status);
      point_cloud_field_t* field = point_cloud_field_new(local_cloud, num_comps);
      size_t size = num_comps * (local_cloud->num_points + local_cloud->num_ghosts);
      MPI_Recv(field->data, (int)size, MPI_REAL_T, 0, rank, comm, &status);
      local_fields[i] = field;
    }
  }

  // Replace the global thingies with the local thingies.
  *cloud = local_cloud;
  for (size_t i = 0; i < num_fields; ++i)
  {
    if (fields[i] != NULL)
      point_cloud_field_free(fields[i]);
    fields[i] = local_fields[i];
  }

  // Clean up.
  if (global_cloud != NULL)
    point_cloud_free(global_cloud);

  STOP_FUNCTION_TIMER();
#endif
}

bool partition_point_cloud(point_cloud_t** cloud,
                           MPI_Comm comm,
                           int* weights,
                           real_t imbalance_tol,
                           point_cloud_field_t** fields,
                           size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((*cloud == NULL) || ((*cloud)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
  {
    STOP_FUNCTION_TIMER();
    return true;
  }

  // Now do the space-filling curve thing.
  int64_t* global_partition;
  if (rank == 0)
  {
    point_cloud_t* cl = *cloud;
    global_partition = partition_points(cl->points, cl->num_points, comm,
                                        weights, imbalance_tol, true);
  }
  else
    global_partition = partition_points(NULL, 0, comm, weights, imbalance_tol, true);
  if (global_partition == NULL) // couldn't balance the load!
  {
    STOP_FUNCTION_TIMER();
    return false;
  }

  // Distribute the point cloud.
  distribute_point_cloud(cloud, comm, global_partition, fields, num_fields);

  // Clean up.
  if (global_partition != NULL)
    polymec_free(global_partition);

  STOP_FUNCTION_TIMER();
  return true;
#else
  return true;
#endif
}

//------------------------------------------------------------------------
//                     Dynamic repartition code below
//------------------------------------------------------------------------

#if POLYMEC_HAVE_MPI

// Fuse a set of subclouds into a single point cloud. Ghost points are not
// permitted. Subclouds are consumed.
static point_cloud_t* fuse_clouds(point_cloud_t** subclouds, size_t num_subclouds)
{
  size_t num_points = 0;
  for (size_t i = 0; i < num_subclouds; ++i)
  {
    ASSERT(subclouds[i]->num_ghosts == 0);
    num_points += subclouds[i]->num_points;
  }

  point_cloud_t* fused_cloud = point_cloud_new(subclouds[0]->comm, num_points);
  int k = 0;
  for (size_t i = 0; i < num_subclouds; ++i)
  {
    for (size_t j = 0; j < subclouds[i]->num_points; ++j, ++k)
      fused_cloud->points[k] = subclouds[i]->points[j];
  }

  // Consume the subclouds.
  for (size_t i = 0; i < num_subclouds; ++i)
  {
    point_cloud_free(subclouds[i]);
    subclouds[i] = NULL;
  }

  return fused_cloud;
}

#endif

bool repartition_point_cloud(point_cloud_t** cloud,
                             int* weights,
                             real_t imbalance_tol,
                             point_cloud_field_t** fields,
                             size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  point_cloud_t* cl = *cloud;
  int nprocs, rank;
  MPI_Comm_size(cl->comm, &nprocs);
  MPI_Comm_rank(cl->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
  {
    STOP_FUNCTION_TIMER();
    return true;
  }

  // Create local partition/load vectors for each process.
  int64_t* local_partition = repartition_points(cl->points, cl->num_points, cl->comm,
                                                weights, imbalance_tol);

  // Redistribute the point cloud and its fields.
  redistribute_point_cloud(cloud, local_partition, fields, num_fields);

  // Clean up.
  polymec_free(local_partition);

  STOP_FUNCTION_TIMER();
  return true;
#else
  return true;
#endif
}

void redistribute_point_cloud(point_cloud_t** cloud,
                              int64_t* local_partition,
                              point_cloud_field_t** fields,
                              size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  point_cloud_t* c = *cloud;

  // Get redistribution data from the local partition vector.
  redistribution_t* redist = redistribution_from_partition((*cloud)->comm, local_partition, (*cloud)->num_points);

  // Post receives for buffer sizes.
  size_t num_receives = redist->receive_procs->size;
  size_t num_sends = redist->send_procs->size;
  int receive_buffer_sizes[num_receives];
  MPI_Request requests[num_receives + num_sends];
  for (size_t p = 0; p < num_receives; ++p)
  {
    int proc = redist->receive_procs->data[p];
    MPI_Irecv(&receive_buffer_sizes[p], 1, MPI_INT, proc, 0, c->comm, &requests[p]);
  }

  // Build point clouds to send to other processes.
  int_unordered_set_t* sent_points = int_unordered_set_new();
  serializer_t* ser = point_cloud_serializer();
  byte_array_t* send_buffers[num_sends];
  for (size_t p = 0; p < num_sends; ++p)
  {
    int proc = redist->send_procs->data[p];
    int_array_t* indices = redist->send_indices[p];
    byte_array_t* bytes = byte_array_new();

    // Add the indices of the cells we are sending.
    for (size_t i = 0; i < indices->size; ++i)
      int_unordered_set_insert(sent_points, indices->data[i]);

    // Create the subcloud to send.
    point_cloud_t* subcloud = create_subcloud(c->comm, c, indices->data, indices->size);

    // Serialize the cloud.
    size_t offset = 0;
    serializer_write(ser, subcloud, bytes, &offset);

    // Pack the relevant contents of the fields into the buffer.
    for (size_t i = 0; i < num_fields; ++i)
    {
      int num_comps = fields[i]->num_components;
      byte_array_write_size_ts(bytes, 1, &subcloud->num_points, &offset);
      real_t field_data[num_comps*subcloud->num_points];
      for (size_t j = 0; j < subcloud->num_points; ++j)
        for (size_t cc = 0; cc < num_comps; ++cc)
          field_data[num_comps*j+cc] = fields[i]->data[num_comps*indices->data[j]+cc];
      byte_array_write_real_ts(bytes, num_comps*subcloud->num_points, field_data, &offset);
    }

    // Send the buffer size.
    MPI_Isend(&bytes->size, 1, MPI_INT, proc, 0, c->comm, &requests[p + num_receives]);

    // Clean up.
    point_cloud_free(subcloud);
    send_buffers[p] = bytes;
  }

  // Wait for the buffer sizes to be transmitted.
  MPI_Status statuses[num_receives + num_sends];
  MPI_Waitall((int)(num_receives + num_sends), requests, statuses);

  // Post receives for the actual messages.
  byte_array_t* receive_buffers[num_receives];
  for (size_t i = 0; i < num_receives; ++i)
  {
    receive_buffers[i] = byte_array_new();
    byte_array_resize(receive_buffers[i], receive_buffer_sizes[i]);
    MPI_Irecv(receive_buffers[i]->data, receive_buffer_sizes[i], MPI_BYTE,
              redist->receive_procs->data[i], 0, c->comm, &requests[i]);
  }

  // Send the actual clouds and wait for receipt.
  for (size_t i = 0; i < num_sends; ++i)
  {
    MPI_Isend(send_buffers[i]->data, (int)send_buffers[i]->size, MPI_BYTE,
              redist->send_procs->data[i], 0, c->comm, &requests[num_receives + i]);
  }
  MPI_Waitall((int)(num_receives + num_sends), requests, statuses);

  // Unpack the clouds.
  size_t receive_offsets[num_receives];
  point_cloud_t* subclouds[1+num_receives];
  for (size_t i = 0; i < num_receives; ++i)
  {
    receive_offsets[i] = 0;
    subclouds[i+1] = serializer_read(ser, receive_buffers[i], &receive_offsets[i]);
  }

  // Clean up the send buffers and the serializer. We still need the
  // receive buffer.
  ser = NULL;
  for (size_t i = 0; i < num_sends; ++i)
    byte_array_free(send_buffers[i]);

  // Construct a local subcloud and store it in subclouds[0]. This subcloud
  // consists of all points not sent to other processes.
  size_t num_local_points = c->num_points - sent_points->size;
  int local_points[num_local_points], j = 0;
  for (size_t i = 0; i < c->num_points; ++i)
  {
    if (!int_unordered_set_contains(sent_points, (int)i))
      local_points[j++] = (int)i;
  }
  subclouds[0] = create_subcloud(c->comm, c, local_points, num_local_points);

  // Fuse all the subclouds into a single point cloud.
  *cloud = fuse_clouds(subclouds, 1+num_receives);

  // Unpack the migrated field data.
  for (size_t i = 0; i < num_fields; ++i)
  {
    int num_comps = fields[i]->num_components;
    point_cloud_field_t* new_field = point_cloud_field_new(*cloud, num_comps);

    // Local portion of the field.
    for (int k = 0; k < num_local_points; ++k)
      for (size_t cc = 0; cc < num_comps; ++cc)
        new_field->data[num_comps*k+cc] = fields[i]->data[num_comps*local_points[k]+cc];
    size_t offset = num_comps*num_local_points;

    // Remote portions.
    for (size_t r = 0; r < num_receives; ++r)
    {
      size_t num_values;
      byte_array_read_size_ts(receive_buffers[r], 1, &num_values, &receive_offsets[r]);
      byte_array_read_real_ts(receive_buffers[r], num_comps*num_values, &(new_field->data[offset]), &receive_offsets[r]);
      offset += num_comps*num_values;
    }

    // Out with the old, in with the new!
    point_cloud_field_free(fields[i]);
    fields[i] = new_field;
  }

  // Clean up the rest.
  for (size_t i = 0; i < num_receives; ++i)
    byte_array_free(receive_buffers[i]);
  int_unordered_set_free(sent_points);
  redistribution_free(redist);
  point_cloud_free(c);
  STOP_FUNCTION_TIMER();
#endif
}

