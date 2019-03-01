// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/blob_exchanger.h"

// A blob buffer is functionally just a big chunk of memory with an associated 
// type. But it also needs to conduct its own message-passing business.
struct blob_buffer_t 
{
  // Sizing factor for blobs.
  int size_factor;

  // Message passing metadata.
  int num_sends;
  int num_receives;
  int num_requests;
  void** send_buffers;
  int* send_buffer_sizes;
  int* dest_procs;
  void** receive_buffers;
  int* receive_buffer_sizes;
  int* source_procs;
  MPI_Request* requests;

  // Blob storage.
  void* storage;
};

void blob_buffer_free(blob_buffer_t* buffer)
{
  if (buffer->dest_procs != NULL)
    polymec_free(buffer->dest_procs);
  if (buffer->send_buffers != NULL)
    polymec_free(buffer->send_buffers);
  if (buffer->send_buffer_sizes != NULL)
    polymec_free(buffer->send_buffer_sizes);
  if (buffer->source_procs != NULL)
    polymec_free(buffer->source_procs);
  if (buffer->receive_buffers != NULL)
    polymec_free(buffer->receive_buffers);
  if (buffer->receive_buffer_sizes != NULL)
    polymec_free(buffer->receive_buffer_sizes);
  if (buffer->requests != NULL)
    polymec_free(buffer->requests);
  polymec_free(buffer);
}

struct blob_exchanger_t 
{
  MPI_Comm comm;

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
  int_array_t** indices_p = exchanger_proc_map_get(map, process);
  int_array_t* indices = NULL;
  if (indices_p != NULL)
    indices = *indices_p;
  else
  {
    indices = int_array_new();
    exchanger_proc_map_insert_with_v_dtor(map, process, indices, int_array_free);
  }
  int_array_append(indices, blob_index);
}

static void blob_exchanger_free(void* context)
{
  blob_exchanger_t* ex = context;
  ptr_array_free(ex->pending_msgs);
}

blob_exchanger_t* blob_exchanger_new(MPI_Comm comm,
                                     blob_exchanger_proc_map_t* send_map,
                                     blob_exchanger_proc_map_t* receive_map,
                                     blob_exchanger_size_map_t* blob_size_map)
{
  blob_exchanger_t* ex = polymec_refcounted_malloc(sizeof(blob_exchanger_t),
                                                   blob_exchanger_free);
  ex->comm = comm;
  return ex;
}

MPI_Comm blob_exchanger_comm(blob_exchanger_t* ex)
{
  return ex->comm;
}

blob_buffer_t* blob_exchanger_create_buffer(blob_exchanger_t* ex,
                                            int size_factor)
{
  ASSERT(size_factor > 0);
  size_t buffer_size = size_factor * ex->buffer_size;
  blob_buffer_t* b = polymec_malloc(sizeof(blob_buffer_t) + buffer_size);
  b->size_factor = size_factor;
  b->storage = b + sizeof(blob_buffer_t);
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
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
  return 0; // FIXME
}

static int blob_exchanger_waitall(blob_exchanger_t* ex, blob_buffer_t* buffer)
{
#if POLYMEC_HAVE_MPI
  // Allocate storage for statuses of sends/receives.
  int num_requests = buffer->num_requests;
  MPI_Status statuses[num_requests];
  
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
    bool expecting_data = (buffer->num_receives > 0);
    bool sent_data = (buffer->num_sends > 0);

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
        int outstanding_send_procs[buffer->num_sends],
            outstanding_send_bytes[buffer->num_sends],
            outstanding_receive_procs[buffer->num_receives],
            outstanding_receive_bytes[buffer->num_receives],
            completed_send_procs[buffer->num_sends],
            completed_send_bytes[buffer->num_sends],
            completed_receive_procs[buffer->num_receives],
            completed_receive_bytes[buffer->num_receives];
        for (int i = 0; i < num_requests; ++i)
        {
          if (!finished[i])
          {
            if (expecting_data && (i < ex->receive_map->size)) // outstanding receive 
            {
              outstanding_receive_procs[num_outstanding_receives] = buffer->source_procs[i];
              outstanding_receive_bytes[num_outstanding_receives] = buffer->receive_buffer_sizes[i] * buffer->data_size;
              ++num_outstanding_receives;
            }
            else if (sent_data) // outstanding send 
            {
              outstanding_send_procs[num_outstanding_sends] = buffer->dest_procs[i - buffer->num_receives];
              outstanding_send_bytes[num_outstanding_sends] = buffer->send_buffer_sizes[i - buffer->num_receives] * buffer->data_size;
              ++num_outstanding_sends;
            }
          }
          else
          {
            if (expecting_data && (i < buffer->num_receives)) // completed receive 
            {
              completed_receive_procs[num_completed_receives] = buffer->source_procs[i];
              completed_receive_bytes[num_completed_receives] = buffer->receive_buffer_sizes[i] * buffer->data_size;
              ++num_completed_receives;
            }
            else if (sent_data) // completed send 
            {
              completed_send_procs[num_completed_sends] = buffer->dest_procs[i - buffer->num_receives];
              completed_send_bytes[num_completed_sends] = buffer->send_buffer_sizes[i - buffer->num_receives] * buffer->data_size;
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
          fprintf(ex->dl_output_stream, "%d: Completed sending data to:\n", ex->rank);
          for (int i = 0; i < num_completed_sends; ++i)
            fprintf(ex->dl_output_stream, "%d:  %d (%d bytes)\n", ex->rank, completed_send_procs[i], completed_send_bytes[i]);
        }
        if (num_completed_receives > 0)
        {
          fprintf(ex->dl_output_stream, "%d: Completed receiving data from:\n", ex->rank);
          for (int i = 0; i < num_completed_receives; ++i)
            fprintf(ex->dl_output_stream, "%d:  %d (%d bytes)\n", ex->rank, completed_receive_procs[i], completed_receive_bytes[i]);
        }
        if (num_outstanding_sends > 0)
        {
          fprintf(ex->dl_output_stream, "%d: Still sending data to:\n", ex->rank);
          for (int i = 0; i < num_outstanding_sends; ++i)
            fprintf(ex->dl_output_stream, "%d:  %d (%d bytes)\n", ex->rank, outstanding_send_procs[i], outstanding_send_bytes[i]);
        }
        if (num_outstanding_receives > 0)
        {
          fprintf(ex->dl_output_stream, "Still expecting data from:\n");
          for (int i = 0; i < num_outstanding_receives; ++i)
            fprintf(ex->dl_output_stream, "  %d (%d bytes)\n", outstanding_receive_procs[i], outstanding_receive_bytes[i]);
        }
        fprintf(ex->dl_output_stream, "%d: Grace period: %g seconds\n", ex->rank, ex->dl_thresh);

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
        if (i < buffer->num_receives)
        {
          // Now we can really get nitty-gritty and try to diagnose the
          // problem carefully! 
          if (statuses[i].MPI_ERROR == MPI_ERR_TRUNCATE)
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n"
                    "(Expected %d bytes)\n", ex->rank, buffer->source_procs[i], statuses[i].MPI_ERROR, 
                    errstr, buffer->receive_buffer_sizes[i]);
          }
          else
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n",
                    ex->rank, buffer->source_procs[i], statuses[i].MPI_ERROR, errstr);
          }
        }
        else 
        {
          fprintf(ex->dl_output_stream, "%d: MPI error sending to %d (%d) %s\n",
                  ex->rank, buffer->dest_procs[i - buffer->num_receives], statuses[i].MPI_ERROR, errstr);
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

void blob_exchanger_finish_exchange(blob_exchanger_t* ex, int token)
{
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);

  // Retrieve the message for the given token.
  blob_buffer_t* buffer = ex->pending_msgs[token];
  exchanger_waitall(ex, buffer);

  // Pull the buffer out of our list of pending messages.
  ex->pending_msgs->data[token] = NULL;
  STOP_FUNCTION_TIMER();
}

bool blob_exchanger_next_send_blob(blob_exchanger_t* ex, 
                                   int* pos, 
                                   int* remote_process, 
                                   int* blob_index, 
                                   size_t* blob_size)
{
  return false;
}

bool blob_exchanger_next_receive_blob(blob_exchanger_t* ex, 
                                      int* pos, 
                                      int* remote_process, 
                                      int* blob_index, 
                                      size_t* blob_size)
{
  return false;
}

void blob_exchanger_copy_in(blob_exchanger_t* ex,
                            int blob_index,
                            int size_factor,
                            void* blob,
                            blob_buffer_t* buffer)
{
  ASSERT(size_factor == buffer->size_factor);
}

void blob_exchanger_copy_out(blob_exchanger_t* ex,
                             blob_buffer_t* buffer,
                             int blob_index,
                             int size_factor,
                             void* blob)
{
  ASSERT(size_factor == buffer->size_factor);
}

bool blob_exchanger_verify(blob_exchanger_t* ex, 
                           void (*handler)(const char* format, ...))
{
  return false;
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
  ex->dl_output_stream = NULL;
}

bool blob_exchanger_deadlock_detection_enabled(blob_exchanger_t* ex)
{
  return (ex->dl_thresh > 0.0);
}

void blob_exchanger_fprintf(blob_exchanger_t* ex, FILE* stream)
{
}

