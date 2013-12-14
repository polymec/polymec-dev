// Copyright (c) 2012-2013, Jeffrey N. Johnson
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


#include "exchanger.h"
#include "core/unordered_map.h"

// This is a record of a single communications channel for an exchanger 
// to send or receive data to or from a remote process.
typedef struct
{
  int num_indices;
  int* indices;
} exchanger_channel_t;

static exchanger_channel_t* exchanger_channel_new(int num_indices, int* indices, bool copy_indices)
{
  ASSERT(num_indices > 0);
  ASSERT(indices != NULL);
  exchanger_channel_t* c = malloc(sizeof(exchanger_channel_t));
  c->num_indices = num_indices;
  if (copy_indices)
  {
    c->indices = malloc(sizeof(int)*num_indices);
    memcpy(c->indices, indices, sizeof(int)*num_indices);
  }
  else
    c->indices = indices; // Assumes control of memory
  return c;
}

static void exchanger_channel_free(exchanger_channel_t* c)
{
  free(c->indices);
  free(c);
}

DEFINE_UNORDERED_MAP(exchanger_map, int, exchanger_channel_t*, int_hash, int_equals)

typedef struct 
{
  MPI_Datatype type;
  int stride;     // Number of data per element
  int data_size;  // Size of a datum.
  int tag;
  int num_sends;
  int num_receives;
  void** send_buffers;
  int* send_buffer_sizes;
  void** receive_buffers;
  int* receive_buffer_sizes;
  MPI_Request* requests;
} mpi_message_t;

#if USE_MPI
static mpi_message_t* mpi_message_new(MPI_Datatype type, int stride, int tag)
{
  ASSERT(stride > 0);
  mpi_message_t* msg = malloc(sizeof(mpi_message_t));
  msg->type = type;
  msg->stride = stride;
  msg->tag = tag;
  if (type == MPI_DOUBLE)
    msg->data_size = sizeof(double);
  else if (type == MPI_FLOAT)
    msg->data_size = sizeof(float);
  else if (type == MPI_INT)
    msg->data_size = sizeof(int);
  else if (type == MPI_LONG)
    msg->data_size = sizeof(long);
  else if (type == MPI_LONG_LONG)
    msg->data_size = sizeof(long long);
  else if (type == MPI_CHAR)
    msg->data_size = sizeof(char);
  return msg;
}

static void mpi_message_pack(mpi_message_t* msg, void* data, 
                             exchanger_map_t* send_map, exchanger_map_t* receive_map)
{
  ASSERT(send_map->size >= 0);
  ASSERT(receive_map->size >= 0);
  int num_sends = send_map->size;
  int num_receives = receive_map->size;
  msg->num_sends = num_sends;
  msg->send_buffer_sizes = malloc(sizeof(int)*msg->num_sends);
  msg->send_buffers = malloc(sizeof(void*)*msg->num_sends);
  msg->num_receives = num_receives;
  msg->receive_buffer_sizes = malloc(sizeof(int)*msg->num_receives);
  msg->receive_buffers = malloc(sizeof(void*)*msg->num_receives);

  int pos = 0, proc, i = 0;
  exchange_channel_t* c;
  while (exchanger_map_next(send_map, &pos, &proc, &c))
  {
    msg->send_buffer_sizes[i] = c->num_indices;
    msg->send_buffers[i] = malloc(c->num_indices*msg->data_size*msg->stride);
    if (msg->type == MPI_DOUBLE)
    {
      double* src = (double*)data;
      double* dest = (double*)(msg->send_buffers[i]);
      for (int j = 0; j < send_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[msg->stride*j+s] = src[send_idx[i][msg->stride*j+s]];
    }
    else if (msg->type == MPI_FLOAT)
    {
      float* src = (float*)data;
      float* dest = (float*)(msg->send_buffers[i]);
      for (int j = 0; j < send_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[msg->stride*j+s] = src[send_idx[i][msg->stride*j+s]];
    }
    else if (msg->type == MPI_INT)
    {
      int* src = (int*)data;
      int* dest = (int*)(msg->send_buffers[i]);
      for (int j = 0; j < send_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[msg->stride*j+s] = src[send_idx[i][msg->stride*j+s]];
    }
    else if (msg->type == MPI_LONG)
    {
      long* src = (long*)data;
      long* dest = (long*)(msg->send_buffers[i]);
      for (int j = 0; j < send_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[msg->stride*j+s] = src[send_idx[i][msg->stride*j+s]];
    }
    else if (msg->type == MPI_LONG_LONG)
    {
      long long* src = (long long*)data;
      long long* dest = (long long*)(msg->send_buffers[i]);
      for (int j = 0; j < send_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[msg->stride*j+s] = src[send_idx[i][msg->stride*j+s]];
    }
    else if (msg->type == MPI_CHAR)
    {
      char* src = (char*)data;
      char* dest = (char*)(msg->send_buffers[i]);
      for (int j = 0; j < send_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[msg->stride*j+s] = src[send_idx[i][msg->stride*j+s]];
    }
    ++i;
  }

  pos = 0, i = 0;
  while (exchanger_map_next(receive_map, &pos, &proc, &c))
  {
    msg->receive_buffer_sizes[i] = c->num_indices;
    msg->receive_buffers[i] = malloc(c->num_indices*msg->data_size*msg->stride);
    msg->requests = malloc((num_sends+num_receives)*sizeof(MPI_Request));
    ++i;
  }
}

static void mpi_message_unpack(mpi_message_t* msg, void* data, 
                               int** receive_idx)
{
  for (int i = 0; i < msg->num_receives; ++i)
  {
    if (msg->type == MPI_DOUBLE)
    {
      double* src = (double*)(msg->send_buffers[i]);
      double* dest = (double*)data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[receive_idx[i][msg->stride*j+s]] = src[msg->stride*j+s];
    }
    else if (msg->type == MPI_FLOAT)
    {
      float* src = (float*)(msg->send_buffers[i]);
      float* dest = (float*)data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[receive_idx[i][msg->stride*j+s]] = src[msg->stride*j+s];
    }
    else if (msg->type == MPI_INT)
    {
      int* src = (int*)(msg->send_buffers[i]);
      int* dest = (int*)data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[receive_idx[i][msg->stride*j+s]] = src[msg->stride*j+s];
    }
    else if (msg->type == MPI_LONG)
    {
      long* src = (long*)(msg->send_buffers[i]);
      long* dest = (long*)data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[receive_idx[i][msg->stride*j+s]] = src[msg->stride*j+s];
    }
    else if (msg->type == MPI_LONG_LONG)
    {
      long long* src = (long long*)(msg->send_buffers[i]);
      long long* dest = (long long*)data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[receive_idx[i][msg->stride*j+s]] = src[msg->stride*j+s];
    }
    else if (msg->type == MPI_CHAR)
    {
      char* src = (char*)(msg->send_buffers[i]);
      char* dest = (char*)data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < msg->stride; ++s)
          dest[receive_idx[i][msg->stride*j+s]] = src[msg->stride*j+s];
    }
  }
}

static void mpi_message_free(mpi_message_t* msg)
{
  for (int i = 0; i < msg->num_sends; ++i)
  {
    if (msg->send_buffers[i] != NULL)
      free(msg->send_buffers[i]);
  }
  free(msg->send_buffers);
  free(msg->send_buffer_sizes);
  for (int i = 0; i < msg->num_receives; ++i)
  {
    if (msg->receive_buffers[i] != NULL)
      free(msg->receive_buffers[i]);
  }
  free(msg->receive_buffers);
  free(msg->receive_buffer_sizes);
  if (msg->requests != NULL)
    free(msg->requests);
  free(msg);
}

static void mpi_message_fprintf(mpi_message_t* msg, FILE* stream)
{
  char typeStr[1024];
  if (msg->type == MPI_DOUBLE)
    strcpy(typeStr, "double");
  else if (msg->type == MPI_FLOAT)
    strcpy(typeStr, "float");
  else if (msg->type == MPI_INT)
    strcpy(typeStr, "int");
  else if (msg->type == MPI_LONG)
    strcpy(typeStr, "long");
  else if (msg->type == MPI_LONG_LONG)
    strcpy(typeStr, "long long");
  else
    strcpy(typeStr, "char");
  fprintf(stream, "Message(datatype = %s, stride = %d, tag = %d):\n",
          typeStr, msg->stride, msg->tag);
  if (msg->num_sends > 0)
    fprintf(stream, "Send buffers:\n");
  for (int i = 0; i < msg->num_sends; ++i)
    fprintf(stream, " %d: %d bytes\n", i, msg->send_buffer_sizes[i]);
  if (msg->num_receives > 0)
    fprintf(stream, "Receive buffers:\n");
  for (int i = 0; i < msg->num_receives; ++i)
    fprintf(stream, " %d: %d bytes\n", i, msg->receive_buffer_sizes[i]);
}
#endif

struct exchanger_t
{
  MPI_Comm comm;
  int rank;

  // Communication maps.
  exchanger_map_t* send_map;
  exchanger_map_t* receive_map;

  int max_send, max_receive;

  // Pending messages.
  int num_pending_msgs;
  int pending_msg_cap;
  mpi_message_t** pending_msgs;
  void** orig_buffers;
  int** transfer_counts;

  // Deadlock detection.
  double dl_thresh;
  int dl_output_rank;
  FILE* dl_output_stream;
};

exchanger_t* exchanger_new(MPI_Comm comm)
{
  exchanger_t* ex = malloc(sizeof(exchanger_t));
  ex->comm = comm;
  MPI_Comm_rank(comm, &(ex->rank));
  ex->dl_thresh = -1.0;
  ex->dl_output_rank = -1;
  ex->dl_output_stream = NULL;
  ex->send_map = exchanger_map_new();
  ex->receive_map = exchanger_map_new();
  ex->pending_msg_cap = 32;
  ex->pending_msgs = calloc(ex->pending_msg_cap, sizeof(mpi_message_t*));
  ex->orig_buffers = calloc(ex->pending_msg_cap, sizeof(void*));
  ex->transfer_counts = calloc(ex->pending_msg_cap, sizeof(int*));

  return ex;
}

static void exchanger_clear(exchanger_t* ex)
{
  exchanger_map_clear(ex->send_map);
  exchanger_map_clear(ex->receive_map);

  free(ex->pending_msgs);
  free(ex->orig_buffers);
  free(ex->transfer_counts);
}

void exchanger_free(exchanger_t* ex)
{
  exchanger_clear(ex);
  exchanger_map_free(ex->send_map);
  exchanger_map_free(ex->receive_map);
  free(ex);
}

static void delete_map_entry(int key, exchanger_channel_t* value)
{
  exchanger_channel_free(value);
}

void exchanger_set_send(exchanger_t* ex, int remote_process, int num_indices, int* indices, bool copy_indices)
{
  exchanger_channel_t* c = exchanger_channel_new(num_indices, indices, copy_indices);
  exchanger_map_insert_with_kv_dtor(ex->send_map, remote_process, c, delete_map_entry);

  if (remote_process > ex->max_send)
    ex->max_send = remote_process;
}

void exchanger_delete_send(exchanger_t* ex, int remote_process)
{
  exchanger_map_delete(ex->send_map, remote_process);

  // Find the maximum rank to which we now send data.
  ex->max_send = -1;
  int pos = 0, proc;
  exchanger_channel_t* c;
  while (exchanger_map_next(ex->send_map, &pos, &proc, &c))
    ex->max_send = (proc > ex->max_send) ? proc : ex->max_send;
}

void exchanger_set_receive(exchanger_t* ex, int remote_process, int num_indices, int* indices, bool copy_indices)
{
  exchanger_channel_t* c = exchanger_channel_new(num_indices, indices, copy_indices);
  exchanger_map_insert_with_kv_dtor(ex->receive_map, remote_process, c, delete_map_entry);

  if (remote_process > ex->max_receive)
    ex->max_receive = remote_process;
}

void exchanger_delete_receive(exchanger_t* ex, int remote_process)
{
  exchanger_map_delete(ex->send_map, remote_process);

  // Find the maximum rank from which we now receive data.
  ex->max_receive = -1;
  int pos = 0, proc;
  exchanger_channel_t* c;
  while (exchanger_map_next(ex->receive_map, &pos, &proc, &c))
    ex->max_receive = (proc > ex->max_receive) ? proc : ex->max_receive;
}

bool exchanger_is_valid(exchanger_t* ex)
{
  // FIXME!!!
  return true; 
}

int exchanger_max_send(exchanger_t* ex)
{
  return ex->max_send;
}

int exchanger_max_receive(exchanger_t* ex)
{
  return ex->max_receive;
}

void exchanger_enable_deadlock_detection(exchanger_t* ex, 
                                         double threshold,
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

void exchanger_disable_deadlock_detection(exchanger_t* ex)
{
  ex->dl_thresh = -1.0;
  ex->dl_output_rank = -1;
  ex->dl_output_stream = NULL;
}

bool exchanger_deadlock_detection_enabled(exchanger_t* ex)
{
  return (ex->dl_thresh > 0.0);
}

void exchanger_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type)
{
  int token = exchanger_start_exchange(ex, data, stride, tag, type);
  exchanger_finish_exchange(ex, token);
}

int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type)
{
#if USE_MPI
  // Create a message for this array.
  mpi_message_t* msg = mpi_message_new(type, stride, tag);
  mpi_message_pack(msg, data, ex->send_map, ex->receive_map);

  // If we are expecting data, post asynchronous receives. 
  for (int i = 0; i < msg->num_receives; ++i)
  {
    int err = MPI_Irecv(msg->receive_buffers[i], 
                        msg->receive_buffer_sizes[i]/msg->data_size, 
                        msg->type, ex->receives[i], msg->tag, ex->comm, 
                        &(msg->requests[i]));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error posting receive from %d: %d\n(%s)\n", 
               ex->rank, ex->receives[i], err, str);
      polymec_error(err_msg);
    }
  }

  // Send the data asynchronously. 
  for (int i = 0; i < msg->num_sends; ++i)
  {
    int err = MPI_Isend(msg->send_buffers[i], 
                        msg->send_buffer_sizes[i]/msg->data_size, 
                        msg->type, ex->sends[i], msg->tag, ex->comm, 
                        &(msg->requests[i + ex->num_receives])); 
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error sending to %d: %d\n(%s)\n", 
               ex->rank, ex->receives[i], err, str);
      polymec_error(err_msg);
    }
  }

  // Allocate a token for the transmission and store the pending message.
  int token = 0;
  while ((token < ex->num_pending_msgs) && (ex->pending_msgs[token] != NULL))
    ++token;
  if (ex->num_pending_msgs == ex->pending_msg_cap)
  {
    ex->pending_msg_cap *= 2;
    ex->pending_msgs = realloc(ex->pending_msgs, ex->pending_msg_cap*sizeof(mpi_message_t));
    ex->orig_buffers = realloc(ex->orig_buffers, ex->pending_msg_cap*sizeof(void*));
    ex->transfer_counts = realloc(ex->transfer_counts, ex->pending_msg_cap*sizeof(int*));
  }
  ex->pending_msgs[token] = msg;
  ex->orig_buffers[token] = data;
  return token;
#else
  return 0;
#endif
}

#if USE_MPI
static int exchanger_waitall(exchanger_t* ex, mpi_message_t* msg)
{
  // Allocate storage for statuses of sends/receives.
  int num_requests = msg->num_sends + msg->num_receives;
  MPI_Status statuses[num_requests];
  
  int err = 0;
  if (ex->dl_thresh <= 0.0)
  {
    // If we're not using deadlock detection, simply call MPI_Waitall. 
    err = MPI_Waitall(num_requests, msg->requests, statuses);
  }
  else
  {
    // Otherwise, we get all fancy.
    int finished[num_requests];
    memset(finished, 0, num_requests*sizeof(int));
    bool expecting_data = (msg->num_receives > 0);
    bool sent_data = (msg->num_sends > 0);

    // Start the deadlock clock.
    double t1 = MPI_Wtime();

    // Now poll the transmissions till they complete.
    bool all_finished;
    do
    {
      all_finished = true;
      for (int i = 0; i < num_requests; ++i)
      {
        if (!finished[i])
        {
          if (MPI_Test(&(msg->requests[i]), &(finished[i]), &(statuses[i])) != MPI_SUCCESS)
            return -1;
          if (!finished[i]) all_finished = false;
        }
      }

      // If the transmissions have finished at this point, we 
      // can break out of the loop. 
      if (all_finished) break;

      // Take a look at the time. 
      double t2 = MPI_Wtime();

      // If we've passed the deadlock threshold, set our error flag and 
      // and gather some diagnostic data. 
      if ((t2 - t1) > ex->dl_thresh)
      {
        // Set a barrier here so that all processes can dump their 
        // data to a comprehensive report. 
        MPI_Barrier(ex->comm);

        // Cancel all communications. 
        for (int i = 0; i < num_requests; ++i)
          MPI_Cancel(&(msg->requests[i]));

        // Now generate a comprehensive report. 
        err = -1;

        int num_outstanding_sends = 0, num_outstanding_receives = 0,
            num_completed_sends = 0, num_completed_receives = 0;
        int outstanding_send_procs[msg->num_sends],
            outstanding_send_bytes[msg->num_sends],
            outstanding_receive_procs[msg->num_receives],
            outstanding_receive_bytes[msg->num_receives],
            completed_send_procs[msg->num_sends],
            completed_send_bytes[msg->num_sends],
            completed_receive_procs[msg->num_receives],
            completed_receive_bytes[msg->num_receives];
        for (int i = 0; i < num_requests; ++i)
        {
          if (!finished[i])
          {
            if (expecting_data && (i < ex->num_receives)) // outstanding receive 
            {
              outstanding_receive_procs[num_outstanding_receives] = ex->receives[i];
              outstanding_receive_bytes[num_outstanding_receives] = msg->receive_buffer_sizes[i];
              ++num_outstanding_receives;
            }
            else if (sent_data) // outstanding send 
            {
              outstanding_send_procs[num_outstanding_sends] = ex->sends[i - msg->num_receives];
              outstanding_send_bytes[num_outstanding_sends] = msg->send_buffer_sizes[i - msg->num_receives];
              ++num_outstanding_sends;
            }
          }
          else
          {
            if (expecting_data && (i < ex->num_receives)) // completed receive 
            {
              completed_receive_procs[num_completed_receives] = ex->receives[i];
              completed_receive_bytes[num_completed_receives] = msg->receive_buffer_sizes[i];
              ++num_completed_receives;
            }
            else if (sent_data) // completed send 
            {
              completed_send_procs[num_completed_sends] = ex->sends[i - msg->num_receives];
              completed_send_bytes[num_completed_sends] = msg->send_buffer_sizes[i - msg->num_receives];
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
          fprintf(ex->dl_output_stream, "Completed sending data to:\n");
          for (int i = 0; i < num_completed_sends; ++i)
            fprintf(ex->dl_output_stream, "  %d (%d bytes)\n", completed_send_procs[i], completed_send_bytes[i]);
        }
        if (num_completed_receives > 0)
        {
          fprintf(ex->dl_output_stream, "Completed receiving data from:\n");
          for (int i = 0; i < num_completed_receives; ++i)
            fprintf(ex->dl_output_stream, "  %d (%d bytes)\n", completed_receive_procs[i], completed_receive_bytes[i]);
        }
        if (num_outstanding_sends > 0)
        {
          fprintf(ex->dl_output_stream, "Still sending data to:\n");
          for (int i = 0; i < num_outstanding_sends; ++i)
            fprintf(ex->dl_output_stream, "  %d (%d bytes)\n", outstanding_send_procs[i], outstanding_send_bytes[i]);
        }
        if (num_outstanding_receives > 0)
        {
          fprintf(ex->dl_output_stream, "Still expecting data from:\n");
          for (int i = 0; i < num_outstanding_receives; ++i)
            fprintf(ex->dl_output_stream, "  %d (%d bytes)\n", outstanding_receive_procs[i], outstanding_receive_bytes[i]);
        }
        fprintf(ex->dl_output_stream, "Grace period: %g seconds\n", ex->dl_thresh);

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
        if (i < ex->num_receives)
        {
          // Now we can really get nitty-gritty and try to diagnose the
          // problem carefully! 
          if (statuses[i].MPI_ERROR == MPI_ERR_TRUNCATE)
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n"
                    "(Expected %d bytes)\n", ex->rank, ex->receives[i], statuses[i].MPI_ERROR, 
                    errstr, msg->receive_buffer_sizes[i]);
          }
          else
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n",
                    ex->rank, ex->receives[i], statuses[i].MPI_ERROR, errstr);
          }
        }
        else 
        {
          fprintf(ex->dl_output_stream, "%d: MPI error sending to %d (%d) %s\n",
                  ex->rank, ex->sends[i - msg->num_receives], statuses[i].MPI_ERROR, errstr);
        }
        return -1;
      }
      // We shouldn't get here. 
    }
  }

  // That's it.
  return err;
}
#endif

void exchanger_finish_exchange(exchanger_t* ex, int token) 
{
#if USE_MPI
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);

  // Retrieve the message for the given token.
  mpi_message_t* msg = ex->pending_msgs[token];
  int err = exchanger_waitall(ex, msg);
  if (err != MPI_SUCCESS) return;

  // Copy the received data into the original array.
  void* orig_buffer = ex->orig_buffers[token];
  mpi_message_unpack(msg, orig_buffer, ex->receive_idx);

  // Pull the message out of our list of pending messages and delete it.
  ex->pending_msgs[token] = NULL;
  ex->orig_buffers[token] = NULL;
  mpi_message_free(msg);
#endif
}

void exchanger_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type)
{
#if USE_MPI
  int token = exchanger_start_transfer(ex, data, count, stride, tag, type);
  exchanger_finish_transfer(ex, token);
#endif
}

int exchanger_start_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type)
{
#if USE_MPI
  // This is the same as an exchange, with an additional datum (the length 
  // of data) stored.
  int token = exchanger_start_exchange(ex, data, stride, tag, type);
  ex->transfer_counts[token] = count; 
  return token;
#else
  return 0;
#endif
}

void exchanger_finish_transfer(exchanger_t* ex, int token)
{
#if USE_MPI
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);

  // This is like exchanger_finish_exchange, but with some rearranging.

  // Retrieve the message for the given token.
  mpi_message_t* msg = ex->pending_msgs[token];
  int* count = ex->transfer_counts[token];
  int err = exchanger_waitall(ex, msg);
  if (err != MPI_SUCCESS) return;

  // Copy the received data into the original array.
  void* orig_buffer = ex->orig_buffers[token];
  mpi_message_unpack(msg, orig_buffer, ex->receive_idx);

  // Now cull the sent elements from the array.
  int num_sent = 0;
  for (int p = 0; p < ex->num_sends; ++p)
  {
    // Move the sent elements to the back.
    int stride = msg->stride;
    if (msg->type == MPI_DOUBLE)
    {
      double* array = (double*)orig_buffer;
      for (int i = 0; i < ex->send_sizes[p]; ++i)
      {
        for (int s = 0; s < stride; ++s)
        {
          array[stride*ex->send_idx[p][i]+s] = array[*count-1-num_sent];
          ++num_sent;
        }
      }
    }
  }

  // Convey the change in count of the transferred data.
  *count -= num_sent;
  ASSERT(*count >= 0);

  // Pull the message out of our list of pending messages and delete it.
  ex->pending_msgs[token] = NULL;
  ex->orig_buffers[token] = NULL;
  ex->transfer_counts[token] = NULL;
  mpi_message_free(msg);
#endif
}

void
exchanger_fprintf(exchanger_t* ex, FILE* stream)
{
#if USE_MPI
  fprintf(stream, "Exchanger(");
  if (ex->comm == MPI_COMM_WORLD)
    fprintf(stream, "MPI_COMM_WORLD");
  else
    fprintf(stream, "[user communicator]");
  fprintf(stream, ", rank %d):", ex->rank);
  if ((ex->num_sends == 0) && (ex->num_receives == 0))
    fprintf(stream, " (no transactions)\n");
  else
  {
    fprintf(stream, "\n");
    for (int i = 0; i < ex->num_sends; ++i)
    {
      fprintf(stream, " [%d] -> [%d: ", ex->rank, ex->sends[i]);
      for (int j = 0; j < ex->send_sizes[i]; ++j)
        fprintf(stream, " %d ", ex->send_idx[i][j]);
      fprintf(stream, "\n");
    }
    for (int i = 0; i < ex->num_receives; ++i)
    {
      fprintf(stream, " [%d] <- [%d: ", ex->rank, ex->receives[i]);
      for (int j = 0; j < ex->receive_sizes[i]; ++j)
        fprintf(stream, " %d ", ex->receive_idx[i][j]);
      fprintf(stream, "\n");
    }
  }
#else
  fprintf(stream, "Exchanger(dummy)\n");
#endif
}


