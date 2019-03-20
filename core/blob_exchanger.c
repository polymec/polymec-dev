// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/array_utils.h"
#include "core/blob_exchanger.h"

// A blob buffer is functionally just a big chunk of memory with an associated
// type. But it also needs to conduct its own message-passing business.
struct blob_buffer_t
{
  // Reference to our blob exchanger.
  blob_exchanger_t* ex;

  // Sizing factor for blobs.
  int size_factor;

  // Message passing metadata.
  int num_requests;
  MPI_Request* requests;

  // Blob storage.
  void* storage;
};

int blob_buffer_size_factor(blob_buffer_t* buffer)
{
  return buffer->size_factor;
}

void blob_buffer_free(blob_buffer_t* buffer)
{
  polymec_free(buffer->storage);
  polymec_free(buffer->requests);
  release_ref(buffer->ex);
  polymec_free(buffer);
}

struct blob_exchanger_t
{
  MPI_Comm comm;
  int rank, nproc;

  // Mapping from processes to indices for blobs sent.
  blob_exchanger_proc_map_t* send_map;
  blob_exchanger_proc_map_t* recv_map;

  // Sizes for blob indices.
  blob_exchanger_size_map_t* blob_sizes;

  // Metadata for creating and manipulating buffers.
  bool does_local_copy;

  int_array_t* send_procs;
  int_array_t* send_proc_offsets;
  int_int_unordered_map_t* send_blob_offsets;
  int max_num_send_blobs;

  int_array_t* recv_procs;
  int_array_t* recv_proc_offsets;
  int_int_unordered_map_t* recv_blob_offsets;
  int max_num_recv_blobs;

  // Base size of buffers created by this blob_exchanger.
  size_t base_buffer_size;

  // Pending messages (buffers).
  ptr_array_t* pending_msgs;

  // Deadlock detection.
  real_t dl_thresh;
  int dl_output_rank;
  FILE* dl_output_stream;
};

void blob_exchanger_proc_map_add_index(blob_exchanger_proc_map_t* map,
                                       int process,
                                       int blob_index)
{
  int_array_t** indices_p = blob_exchanger_proc_map_get(map, process);
  int_array_t* indices = NULL;
  if (indices_p != NULL)
    indices = *indices_p;
  else
  {
    indices = int_array_new();
    blob_exchanger_proc_map_insert_with_v_dtor(map, process, indices,
                                               int_array_free);
  }
  int_array_append(indices, blob_index);
}

static void blob_exchanger_free(void* context)
{
  blob_exchanger_t* ex = context;
  ptr_array_free(ex->pending_msgs);

  int_array_free(ex->send_procs);
  int_array_free(ex->send_proc_offsets);
  int_int_unordered_map_free(ex->send_blob_offsets);

  int_array_free(ex->recv_procs);
  int_array_free(ex->recv_proc_offsets);
  int_int_unordered_map_free(ex->recv_blob_offsets);

  blob_exchanger_proc_map_free(ex->send_map);
  blob_exchanger_proc_map_free(ex->recv_map);
  blob_exchanger_size_map_free(ex->blob_sizes);
}

static void compute_offsets(blob_exchanger_t* ex)
{
  int_array_append(ex->send_proc_offsets, 0);
  ex->max_num_send_blobs = 0;
  int pos = 0, proc;
  int_array_t* indices;
  while (blob_exchanger_proc_map_next(ex->send_map, &pos, &proc, &indices))
  {
    if (proc == ex->rank)
      ex->does_local_copy = true;
    int_array_append(ex->send_procs, proc);
    size_t p = ex->send_proc_offsets->size-1;
    int_array_append(ex->send_proc_offsets, ex->send_proc_offsets->data[p]);
    ex->max_num_send_blobs = MAX(ex->max_num_send_blobs, (int)indices->size);
    for (size_t i = 0; i < indices->size; ++i)
    {
      int blob_index = indices->data[i];
      int_int_unordered_map_insert(ex->send_blob_offsets, blob_index,
                                   ex->send_proc_offsets->data[p+1]);
      size_t blob_size = *blob_exchanger_size_map_get(ex->blob_sizes, blob_index);
      ex->send_proc_offsets->data[p+1] += blob_size;
    }
  }

  // Receive buffers start right after send buffers.
  if (ex->send_proc_offsets->size > 0)
  {
    int_array_append(ex->recv_proc_offsets,
                     ex->send_proc_offsets->data[ex->send_proc_offsets->size-1]);
  }
  else
    int_array_append(ex->recv_proc_offsets, 0);

  ex->max_num_recv_blobs = 0;
  pos = 0;
  while (blob_exchanger_proc_map_next(ex->recv_map, &pos, &proc, &indices))
  {
    int_array_append(ex->recv_procs, proc);
    size_t p = ex->recv_proc_offsets->size-1;
    int_array_append(ex->recv_proc_offsets, ex->recv_proc_offsets->data[p]);
    ex->max_num_recv_blobs = MAX(ex->max_num_recv_blobs, (int)indices->size);
    for (size_t i = 0; i < indices->size; ++i)
    {
      int blob_index = indices->data[i];
      int_int_unordered_map_insert(ex->recv_blob_offsets, blob_index,
                                   ex->recv_proc_offsets->data[p+1]);
      size_t blob_size = *blob_exchanger_size_map_get(ex->blob_sizes, blob_index);
      ex->recv_proc_offsets->data[p+1] += blob_size;
    }
  }

  ex->base_buffer_size = ex->recv_proc_offsets->data[ex->recv_procs->size];
}

blob_exchanger_t* blob_exchanger_new(MPI_Comm comm,
                                     blob_exchanger_proc_map_t* send_map,
                                     blob_exchanger_proc_map_t* receive_map,
                                     blob_exchanger_size_map_t* blob_size_map)
{
  int rank, nproc;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nproc);

  // Validate the maps against one another. There must be a blob size for every
  // index found in the send and receive maps.
#ifndef NDEBUG
  size_t local_bytes_sent = 0;
  int pos = 0, proc;
  int_array_t* indices;
  while (blob_exchanger_proc_map_next(send_map, &pos, &proc, &indices))
  {
    for (size_t i = 0; i < indices->size; ++i)
    {
      int index = indices->data[i];
      ASSERT(blob_exchanger_size_map_contains(blob_size_map, index));
      ASSERT(proc >= 0);
      ASSERT(proc < nproc);
      if (proc == rank)
        local_bytes_sent += *blob_exchanger_size_map_get(blob_size_map, index);
    }
  }

  size_t local_bytes_received = 0;
  pos = 0;
  while (blob_exchanger_proc_map_next(receive_map, &pos, &proc, &indices))
  {
    for (size_t i = 0; i < indices->size; ++i)
    {
      int index = indices->data[i];
      ASSERT(blob_exchanger_size_map_contains(blob_size_map, index));
      ASSERT(proc >= 0);
      if (proc == rank)
        local_bytes_received += *blob_exchanger_size_map_get(blob_size_map, index);
    }
  }

  // Make sure any local copies are consistent.
  ASSERT(local_bytes_sent == local_bytes_received);
#endif

  blob_exchanger_t* ex = polymec_refcounted_malloc(sizeof(blob_exchanger_t),
                                                   blob_exchanger_free);
  ex->comm = comm;
  ex->rank = rank;
  ex->nproc = nproc;
  ex->send_map = send_map;
  ex->recv_map = receive_map;
  ex->blob_sizes = blob_size_map;
  ex->send_procs = int_array_new();
  ex->send_proc_offsets = int_array_new();
  ex->send_blob_offsets = int_int_unordered_map_new();
  ex->recv_procs = int_array_new();
  ex->recv_proc_offsets = int_array_new();
  ex->recv_blob_offsets = int_int_unordered_map_new();
  ex->pending_msgs = ptr_array_new();
  ex->does_local_copy = false;
  ex->dl_thresh = -1.0;
  ex->dl_output_rank = -1;
  ex->dl_output_stream = NULL;

  // Compute offsets for the send/receive blobs.
  compute_offsets(ex);

  return ex;
}

MPI_Comm blob_exchanger_comm(blob_exchanger_t* ex)
{
  return ex->comm;
}

size_t blob_exchanger_blob_size(blob_exchanger_t* ex, int blob_index)
{
  size_t* size_p = blob_exchanger_size_map_get(ex->blob_sizes, blob_index);
  if (size_p != NULL)
    return *size_p;
  else
    return 0;
}

blob_buffer_t* blob_exchanger_create_buffer(blob_exchanger_t* ex,
                                            int size_factor)
{
  ASSERT(size_factor > 0);
  blob_buffer_t* b = polymec_malloc(sizeof(blob_buffer_t));
  retain_ref(ex);
  b->ex = ex;
  b->size_factor = size_factor;
  b->num_requests = (int)((ex->does_local_copy) ?
                           ex->send_procs->size + ex->recv_procs->size - 2 :
                           ex->send_procs->size + ex->recv_procs->size);
  b->requests = polymec_malloc(sizeof(MPI_Request) * b->num_requests);
  size_t buffer_size = size_factor * ex->base_buffer_size;
  b->storage = polymec_malloc(sizeof(char)*buffer_size);
  return b;
}

void blob_exchanger_exchange(blob_exchanger_t* ex,
                             int tag,
                             blob_buffer_t* buffer)
{
  START_FUNCTION_TIMER();
  int token = blob_exchanger_start_exchange(ex, tag, buffer);
  blob_exchanger_finish_exchange(ex, token);
  STOP_FUNCTION_TIMER();
}

int blob_exchanger_start_exchange(blob_exchanger_t* ex,
                                  int tag,
                                  blob_buffer_t* buffer)
{
  ASSERT(buffer->ex == ex);
  START_FUNCTION_TIMER();

  // Append this buffer to the end of our list of pending transactions.
  int token = 0;
  while ((token < ex->pending_msgs->size) &&
         (ex->pending_msgs->data[token] != NULL))
    ++token;
  if (token == (int)ex->pending_msgs->size)
    ptr_array_append(ex->pending_msgs, buffer);
  else
    ex->pending_msgs->data[token] = buffer;

  size_t size_factor = buffer->size_factor;

  // Post receives for the messages.
#if POLYMEC_HAVE_MPI
  int r = 0;
#endif
  for (size_t p = 0; p < ex->recv_procs->size; ++p)
  {
    int proc = ex->recv_procs->data[p];
    size_t r_offset = size_factor * ex->recv_proc_offsets->data[p];
    char* data = &(((char*)buffer->storage)[r_offset]);
    if (proc != ex->rank)
    {
#if POLYMEC_HAVE_MPI
      size_t size = size_factor *
                    (ex->recv_proc_offsets->data[p+1] -
                     ex->recv_proc_offsets->data[p]);
      int err = MPI_Irecv(data, (int)size, MPI_BYTE, proc, tag, ex->comm,
                          &(buffer->requests[r]));
      if (err != MPI_SUCCESS)
      {
        int resultlen;
        char str[MPI_MAX_ERROR_STRING];
        MPI_Error_string(err, str, &resultlen);
        char err_msg[1024];
        snprintf(err_msg, 1024, "%d: MPI Error posting receive from %d: %d\n(%s)\n",
                 ex->rank, proc, err, str);
        polymec_error(err_msg);
      }
      ++r;
#endif
    }
    else // local copy
    {
      size_t p1 = int_lsearch(ex->send_proc_offsets->data,
                              ex->send_proc_offsets->size, ex->rank) -
                  ex->send_proc_offsets->data;
      size_t send_size = size_factor *
                         (ex->send_proc_offsets->data[p1+1] -
                          ex->send_proc_offsets->data[p1]);
      size_t s_offset = size_factor * ex->send_proc_offsets->data[p1];
      char* send_data = &(((char*)buffer->storage)[s_offset]);
      memmove(data, send_data, send_size);
    }
  }

#if POLYMEC_HAVE_MPI
  // Schedule the messages for sending.
  for (size_t p = 0; p < ex->send_procs->size; ++p)
  {
    int proc = ex->send_procs->data[p];
    size_t size = size_factor *
                  (ex->send_proc_offsets->data[p+1] -
                   ex->send_proc_offsets->data[p]);
    size_t s_offset = size_factor * ex->send_proc_offsets->data[p];
    char* data = &(((char*)buffer->storage)[s_offset]);
    if (proc != ex->rank)
    {
      int err = MPI_Isend(data, (int)size, MPI_BYTE, proc, tag, ex->comm,
                          &buffer->requests[r]);
      if (err != MPI_SUCCESS)
      {
        int resultlen;
        char str[MPI_MAX_ERROR_STRING];
        MPI_Error_string(err, str, &resultlen);
        char err_msg[1024];
        snprintf(err_msg, 1024, "%d: MPI Error sending from %d: %d\n(%s)\n",
                 ex->rank, proc, err, str);
        polymec_error(err_msg);
      }
      ++r;
    }
  }
  ASSERT(r == buffer->num_requests);
#endif

  STOP_FUNCTION_TIMER();
  return token;
}

static int blob_exchanger_waitall(blob_exchanger_t* ex, blob_buffer_t* buffer)
{
#if POLYMEC_HAVE_MPI
  // Allocate storage for statuses of sends/receives.
  int num_requests = buffer->num_requests;
  MPI_Status statuses[num_requests];

  int num_sends = (int)ex->send_procs->size;
  int num_receives = (int)ex->recv_procs->size;
  if (ex->does_local_copy)
  {
    --num_sends;
    --num_receives;
  }
  ASSERT(num_requests == num_sends + num_receives);

  // If we're not actually doing remote communication, there's nothing
  // to do here.
  if (num_requests == 0)
    return 0;

  int err = 0;
  if (ex->dl_thresh <= 0.0)
  {
    // If we're not using deadlock detection, simply call MPI_Waitall.
    err = MPI_Waitall(num_requests, buffer->requests, statuses);
  }
  else
  {
    // Otherwise, we get all fancy.
    int finished[num_requests];
    memset(finished, 0, num_requests*sizeof(int));
    bool expecting_data = (num_receives > 0);
    bool sent_data = (num_sends > 0);

    // Start the deadlock clock.
    real_t t1 = (real_t)MPI_Wtime();

    // Now poll the transmissions till they complete.
    bool all_finished;
    do
    {
      all_finished = true;
      for (int i = 0; i < num_requests; ++i)
      {
        if (!finished[i])
        {
          if (MPI_Test(&(buffer->requests[i]), &(finished[i]), &(statuses[i])) != MPI_SUCCESS)
            return -1;
          if (!finished[i]) all_finished = false;
        }
      }

      // If the transmissions have finished at this point, we
      // can break out of the loop.
      if (all_finished) break;

      // Take a look at the time.
      real_t t2 = (real_t)MPI_Wtime();

      // If we've passed the deadlock threshold, set our error flag and
      // and gather some diagnostic data.
      if ((t2 - t1) > ex->dl_thresh)
      {
        // Cancel all unfinished communications.
        for (int i = 0; i < num_requests; ++i)
        {
          if (!finished[i])
            MPI_Cancel(&(buffer->requests[i]));
        }

        // Now generate a comprehensive report.
        err = -1;

        int num_outstanding_sends = 0, num_outstanding_receives = 0,
            num_completed_sends = 0, num_completed_receives = 0;
        int outstanding_send_procs[num_sends],
            outstanding_send_bytes[num_sends],
            outstanding_recv_procs[num_receives],
            outstanding_recv_bytes[num_receives],
            completed_send_procs[num_sends],
            completed_send_bytes[num_sends],
            completed_recv_procs[num_receives],
            completed_recv_bytes[num_receives];
        for (int r = 0; r < num_requests; ++r)
        {
          if (!finished[r])
          {
            if (expecting_data && (r < ex->recv_map->size)) // outstanding receive
            {
              outstanding_recv_procs[num_outstanding_receives] = ex->recv_procs->data[r];
              int size = ex->recv_proc_offsets->data[r+1] - ex->recv_proc_offsets->data[r];
              outstanding_recv_bytes[num_outstanding_receives] = size;
              ++num_outstanding_receives;
            }
            else if (sent_data) // outstanding send
            {
              int s = r - (int)ex->recv_procs->size;
              outstanding_send_procs[num_outstanding_sends] = ex->send_procs->data[s];
              int size = ex->send_proc_offsets->data[s+1] - ex->send_proc_offsets->data[s];
              outstanding_send_bytes[num_outstanding_sends] = size;
              ++num_outstanding_sends;
            }
          }
          else
          {
            if (expecting_data && (r < num_receives)) // completed receive
            {
              completed_recv_procs[num_completed_receives] = ex->recv_procs->data[r];
              int size = ex->recv_proc_offsets->data[r+1] - ex->recv_proc_offsets->data[r];
              completed_recv_bytes[num_completed_receives] = size;
              ++num_completed_receives;
            }
            else if (sent_data) // completed send
            {
              int s = r - (int)ex->recv_procs->size;
              completed_send_procs[num_completed_sends] = ex->send_procs->data[s];
              int size = ex->send_proc_offsets->data[s+1] - ex->send_proc_offsets->data[s];
              completed_send_bytes[num_completed_sends] = size;
              ++num_completed_sends;
            }
          }
        }

        // At this point, there must be at least one uncompleted
        // send and/or receive.
        ASSERT((num_outstanding_sends > 0) || (num_outstanding_receives > 0));

        // Format the report.
        fprintf(ex->dl_output_stream, "%d: MPI Deadlock:\n", ex->rank);
        if (num_completed_sends > 0)
        {
          fprintf(ex->dl_output_stream, "%d: Completed sending data to:\n",
                  ex->rank);
          for (int i = 0; i < num_completed_sends; ++i)
          {
            fprintf(ex->dl_output_stream, "%d:  %d (%d bytes)\n", ex->rank,
                    completed_send_procs[i], completed_send_bytes[i]);
          }
        }
        if (num_completed_receives > 0)
        {
          fprintf(ex->dl_output_stream, "%d: Completed receiving data from:\n",
                  ex->rank);
          for (int i = 0; i < num_completed_receives; ++i)
          {
            fprintf(ex->dl_output_stream, "%d:  %d (%d bytes)\n", ex->rank,
                    completed_recv_procs[i], completed_recv_bytes[i]);
          }
        }
        if (num_outstanding_sends > 0)
        {
          fprintf(ex->dl_output_stream, "%d: Still sending data to:\n",
                  ex->rank);
          for (int i = 0; i < num_outstanding_sends; ++i)
          {
            fprintf(ex->dl_output_stream, "%d:  %d (%d bytes)\n", ex->rank,
                    outstanding_send_procs[i], outstanding_send_bytes[i]);
          }
        }
        if (num_outstanding_receives > 0)
        {
          fprintf(ex->dl_output_stream, "Still expecting data from:\n");
          for (int i = 0; i < num_outstanding_receives; ++i)
          {
            fprintf(ex->dl_output_stream, "  %d (%d bytes)\n",
                    outstanding_recv_procs[i], outstanding_recv_bytes[i]);
          }
        }
        fprintf(ex->dl_output_stream, "%d: Grace period: %g seconds\n",
                ex->rank, ex->dl_thresh);

        // Bug out.
        return -1;
      }
      // Otherwise, slog onward.
    }
    while (!all_finished && (err == 0));
  }

  // If the status buffer contains any errors, check it out.
  if (err == MPI_ERR_IN_STATUS)
  {
    char errstr[MPI_MAX_ERROR_STRING];
    int errlen;
    for (int i = 0; i < num_requests; ++i)
    {
      if (statuses[i].MPI_ERROR != MPI_SUCCESS)
      {
        MPI_Error_string(statuses[i].MPI_ERROR, errstr, &errlen);
        if (i < num_receives)
        {
          // Now we can really get nitty-gritty and try to diagnose the
          // problem carefully!
          if (statuses[i].MPI_ERROR == MPI_ERR_TRUNCATE)
          {
            int size = ex->recv_proc_offsets->data[i+1] -
                       ex->recv_proc_offsets->data[i];
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n"
                    "(Expected %d bytes)\n", ex->rank, statuses[i].MPI_SOURCE,
                    statuses[i].MPI_ERROR, errstr, size);
          }
          else
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n",
                    ex->rank, statuses[i].MPI_SOURCE, statuses[i].MPI_ERROR,
                    errstr);
          }
        }
        else
        {
          fprintf(ex->dl_output_stream, "%d: MPI error sending to %d (%d) %s\n",
                  ex->rank, ex->send_procs->data[i - num_receives],
                  statuses[i].MPI_ERROR, errstr);
        }
        return -1;
      }
      // We shouldn't get here.
    }
  }

  // That's it.
  return err;
#else
  return 0;
#endif
}

bool blob_exchanger_finish_exchange(blob_exchanger_t* ex, int token)
{
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  if (token >= (int)ex->pending_msgs->size)
    return false;

  // Retrieve the message for the given token.
  blob_buffer_t* buffer = ex->pending_msgs->data[token];
  ASSERT(buffer != NULL);
  blob_exchanger_waitall(ex, buffer);

  // Pull the buffer out of our list of pending messages.
  ex->pending_msgs->data[token] = NULL;
  STOP_FUNCTION_TIMER();
  return true;
}

bool blob_exchanger_next_send_blob(blob_exchanger_t* ex,
                                   int* pos,
                                   int* remote_process,
                                   int* blob_index,
                                   size_t* blob_size)
{
  int pos1 = *pos / ex->max_num_send_blobs;
  int_array_t* indices;
  bool result = blob_exchanger_proc_map_next(ex->send_map, &pos1,
                                             remote_process, &indices);
  if (!result)
    return false;

  int i = *pos % ex->max_num_send_blobs;
  *blob_index = indices->data[i];
  *blob_size = *blob_exchanger_size_map_get(ex->blob_sizes, *blob_index);
  if (i < (int)indices->size)
    ++(*pos);
  else
    *pos = pos1 * ex->max_num_send_blobs;
  return true;
}

bool blob_exchanger_next_receive_blob(blob_exchanger_t* ex,
                                      int* pos,
                                      int* remote_process,
                                      int* blob_index,
                                      size_t* blob_size)
{
  int pos1 = *pos / ex->max_num_recv_blobs;
  int_array_t* indices;
  bool result = blob_exchanger_proc_map_next(ex->recv_map, &pos1,
                                             remote_process, &indices);
  if (!result)
    return false;

  int i = *pos % ex->max_num_recv_blobs;
  *blob_index = indices->data[i];
  *blob_size = *blob_exchanger_size_map_get(ex->blob_sizes, *blob_index);
  if (i < (int)indices->size)
    ++(*pos);
  else
    *pos = pos1 * ex->max_num_recv_blobs;
  return true;
}

bool blob_exchanger_copy_in(blob_exchanger_t* ex,
                            int blob_index,
                            void* blob,
                            blob_buffer_t* buffer)
{
  ASSERT(buffer->ex == ex);
  size_t size_factor = buffer->size_factor;
  int* offset_p = int_int_unordered_map_get(ex->send_blob_offsets, blob_index);
  if (offset_p != NULL)
  {
    size_t offset = size_factor * (*offset_p);
    size_t size = *blob_exchanger_size_map_get(ex->blob_sizes, blob_index);
    memcpy(&(((char*)buffer->storage)[offset]), blob, size*size_factor);
    return true;
  }
  else
    return false;
}

bool blob_exchanger_copy_out(blob_exchanger_t* ex,
                             blob_buffer_t* buffer,
                             int blob_index,
                             void* blob)
{
  ASSERT(buffer->ex == ex);
  size_t size_factor = buffer->size_factor;
  int* offset_p = int_int_unordered_map_get(ex->recv_blob_offsets, blob_index);
  if (offset_p != NULL)
  {
    size_t offset = size_factor * (*offset_p);
    size_t size = *blob_exchanger_size_map_get(ex->blob_sizes, blob_index);
    memcpy(blob, &(((char*)buffer->storage)[offset]), size*size_factor);
    return true;
  }
  else
    return false;
}

bool blob_exchanger_is_valid(blob_exchanger_t* ex, char** reason)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  log_debug("blob_exchanger_verify: Checking connectivity.");

  static char _reason[1025];

  // An exchanger is valid/consistent iff the number of blobs that
  // are exchanged between any two processors are agreed upon between those
  // two processors.

  // First question: do our neighbors agree with us about being our
  // neighbors?

  // Tally up the number of indices we're sending to and receiving from
  // every other process in the communicator.
  int pos = 0, proc;
  int_array_t* indices;
  int my_neighbors[2*ex->nproc];
  memset(my_neighbors, 0, 2 * sizeof(int) * ex->nproc);
  while (blob_exchanger_proc_map_next(ex->send_map, &pos, &proc, &indices))
    my_neighbors[2*proc] += (int)indices->size;
  pos = 0;
  while (blob_exchanger_proc_map_next(ex->recv_map, &pos, &proc, &indices))
    my_neighbors[2*proc+1] += (int)indices->size;

  // Do an all-to-all exchange to get everyone's votes on who is whose
  // neighbor.
  int* neighbors_for_proc = polymec_malloc(sizeof(int)*2*ex->nproc*ex->nproc);
  MPI_Allgather(my_neighbors, 2*ex->nproc, MPI_INT,
                neighbors_for_proc, 2*ex->nproc, MPI_INT, ex->comm);

  for (int p = 0; p < ex->nproc; ++p)
  {
    int num_im_sending = my_neighbors[2*p];
    int num_theyre_receiving = neighbors_for_proc[2*(p*ex->nproc+ex->rank)+1];
    if (num_im_sending != num_theyre_receiving)
    {
      polymec_free(neighbors_for_proc);
      snprintf(_reason, 1024, "Proc %d is sending %d blobs to proc %d,"
               " which is expecting %d blobs.",
               ex->rank, num_im_sending, p, num_theyre_receiving);
      if (reason != NULL)
        *reason = _reason;
      STOP_FUNCTION_TIMER();
      return false;
    }

    int num_im_receiving = my_neighbors[2*p+1];
    int num_theyre_sending = neighbors_for_proc[2*(p*ex->nproc+ex->rank)];
    if (num_im_receiving != num_theyre_sending)
    {
      polymec_free(neighbors_for_proc);
      snprintf(_reason, 1024, "Proc %d is sending %d blobs to proc %d,"
              " which is expecting %d blobs.", p, num_theyre_sending,
              ex->rank, num_im_receiving);
      if (reason != NULL)
        *reason = _reason;
      STOP_FUNCTION_TIMER();
      return false;
    }
  }
  polymec_free(neighbors_for_proc);
  log_debug("blob_exchanger_verify: Connectivity verified successfully.");
  STOP_FUNCTION_TIMER();
#endif
  return true;
}

void blob_exchanger_enable_deadlock_detection(blob_exchanger_t* ex,
                                              real_t threshold,
                                              int output_rank,
                                              FILE* stream)
{
  ASSERT(threshold > 0.0);
  ASSERT(output_rank >= 0);
  ASSERT(stream != NULL);
  ex->dl_thresh = threshold;
  ex->dl_output_rank = output_rank;
  ex->dl_output_stream = stream;
}

void blob_exchanger_disable_deadlock_detection(blob_exchanger_t* ex)
{
  ex->dl_thresh = -1.0;
  ex->dl_output_rank = -1;
  if ((ex->dl_output_stream != NULL) &&
      (ex->dl_output_stream != stdout) &&
      (ex->dl_output_stream != stderr))
    fclose(ex->dl_output_stream);
  ex->dl_output_stream = NULL;
}

bool blob_exchanger_deadlock_detection_enabled(blob_exchanger_t* ex)
{
  return (ex->dl_thresh > 0.0);
}

void blob_exchanger_fprintf(blob_exchanger_t* ex, FILE* stream)
{
  fprintf(stream, "blob_exchanger(");
  if (ex->comm == MPI_COMM_WORLD)
    fprintf(stream, "MPI_COMM_WORLD");
  else if (ex->comm == MPI_COMM_SELF)
    fprintf(stream, "MPI_COMM_SELF");
  else
    fprintf(stream, "[user communicator]");
  fprintf(stream, ", rank %d):", ex->rank);
  if ((ex->send_map->size == 0) && (ex->recv_map->size == 0))
    fprintf(stream, " (no transactions)\n");
  else
  {
    fprintf(stream, "\n");
    int pos = 0, proc;
    int_array_t* indices;
    while (blob_exchanger_proc_map_next(ex->send_map, &pos, &proc, &indices))
    {
      fprintf(stream, " [%d] -> [%d]: ", ex->rank, proc);
      for (size_t j = 0; j < indices->size; ++j)
      {
        int index = indices->data[j];
        fprintf(stream, " %d(%dB)", index,
                (int)(*blob_exchanger_size_map_get(ex->blob_sizes, index)));
      }
      fprintf(stream, " ... (%d blobs)\n", (int)indices->size);
    }
    pos = 0;
    while (blob_exchanger_proc_map_next(ex->recv_map, &pos, &proc, &indices))
    {
      fprintf(stream, " [%d] <- [%d]: ", ex->rank, proc);
      for (size_t j = 0; j < indices->size; ++j)
      {
        int index = indices->data[j];
        fprintf(stream, " %d(%dB)", index,
                (int)(*blob_exchanger_size_map_get(ex->blob_sizes, index)));
      }
      fprintf(stream, " ... (%d blobs)\n", (int)indices->size);
    }
  }
}

