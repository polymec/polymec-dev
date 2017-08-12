// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "exchanger.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/array.h"
#include "core/timer.h"

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
  exchanger_channel_t* c = polymec_malloc(sizeof(exchanger_channel_t));
  c->num_indices = num_indices;
  if (copy_indices)
  {
    c->indices = polymec_malloc(sizeof(int)*num_indices);
    memcpy(c->indices, indices, sizeof(int)*num_indices);
  }
  else
    c->indices = indices; // Assumes control of memory
  return c;
}

static void exchanger_channel_free(exchanger_channel_t* c)
{
  polymec_free(c->indices);
  polymec_free(c);
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
  int num_requests;
  void** send_buffers;
  int* send_buffer_sizes;
  int* dest_procs;
  void** receive_buffers;
  int* receive_buffer_sizes;
  int* source_procs;
  MPI_Request* requests;
} mpi_message_t;

#if POLYMEC_HAVE_MPI
static size_t mpi_size(MPI_Datatype type)
{
  size_t size = 0;
  if (type == MPI_REAL_T)
    size = sizeof(real_t);
  else if (type == MPI_DOUBLE)
    size = sizeof(double);
  else if (type == MPI_FLOAT)
    size = sizeof(float);
  else if (type == MPI_INT)
    size = sizeof(int);
  else if (type == MPI_LONG)
    size = sizeof(long);
  else if (type == MPI_LONG_LONG)
    size = sizeof(long long);
  else if (type == MPI_UINT64_T)
    size = sizeof(uint64_t);
  else if (type == MPI_INT64_T)
    size = sizeof(int64_t);
  else if (type == MPI_CHAR)
    size = sizeof(char);
  else if (type == MPI_BYTE)
    size = sizeof(uint8_t);
  else 
    polymec_error("Unsupported MPI data type used.");
  return size;
}

static mpi_message_t* mpi_message_new(MPI_Datatype type, int stride, int tag)
{
  ASSERT(stride > 0);
  mpi_message_t* msg = polymec_malloc(sizeof(mpi_message_t));
  msg->type = type;
  msg->stride = stride;
  msg->tag = tag;
  msg->data_size = (int)mpi_size(type);
  msg->num_sends = 0;
  msg->num_receives = 0;
  msg->num_requests = 0;
  msg->requests = NULL;
  msg->send_buffers = NULL;
  msg->send_buffer_sizes = NULL;
  msg->dest_procs = NULL;
  msg->receive_buffers = NULL;
  msg->receive_buffer_sizes = NULL;
  msg->source_procs = NULL;
  return msg;
}

#define HANDLE_PACKING(msg, mpi_type, c_type) \
  if (msg->type == mpi_type) \
  { \
    c_type* src = data; \
    c_type* dest = msg->send_buffers[i]; \
    for (int j = 0; j < msg->send_buffer_sizes[i]; ++j) \
      for (int s = 0; s < stride; ++s) \
        dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s]; \
  }

static void mpi_message_pack(mpi_message_t* msg, 
                             void* data, 
                             ssize_t send_offset,
                             exchanger_map_t* send_map, 
                             exchanger_map_t* receive_map)
{
  ASSERT(send_map->size >= 0);
  ASSERT(receive_map->size >= 0);
  int num_sends = send_map->size;
  int num_receives = receive_map->size;
  msg->num_sends = num_sends;
  msg->dest_procs = polymec_malloc(sizeof(int)*msg->num_sends);
  msg->send_buffer_sizes = polymec_malloc(sizeof(int)*msg->num_sends);
  msg->send_buffers = polymec_malloc(sizeof(void*)*msg->num_sends);
  msg->num_receives = num_receives;
  msg->source_procs = polymec_malloc(sizeof(int)*msg->num_receives);
  msg->receive_buffer_sizes = polymec_malloc(sizeof(int)*msg->num_receives);
  msg->receive_buffers = polymec_malloc(sizeof(void*)*msg->num_receives);

  int pos = 0, proc, i = 0;
  exchanger_channel_t* c;
  while (exchanger_map_next(send_map, &pos, &proc, &c))
  {
    int num_send_indices = c->num_indices;
    int* send_indices = c->indices;
    int stride = msg->stride;
    msg->dest_procs[i] = proc;
    msg->send_buffer_sizes[i] = num_send_indices;
    msg->send_buffers[i] = polymec_malloc(num_send_indices*msg->data_size*stride);

    HANDLE_PACKING(msg, MPI_REAL_T, real_t)
    HANDLE_PACKING(msg, MPI_FLOAT, float)
    HANDLE_PACKING(msg, MPI_INT, int)
    HANDLE_PACKING(msg, MPI_LONG, long)
    HANDLE_PACKING(msg, MPI_LONG_LONG, long long)
    HANDLE_PACKING(msg, MPI_UINT64_T, uint64_t)
    HANDLE_PACKING(msg, MPI_INT64_T, int64_t)
    HANDLE_PACKING(msg, MPI_CHAR, char)
    HANDLE_PACKING(msg, MPI_BYTE, uint8_t)
    ++i;
  }

  pos = 0; i = 0;
  while (exchanger_map_next(receive_map, &pos, &proc, &c))
  {
    msg->receive_buffer_sizes[i] = c->num_indices;
    msg->receive_buffers[i] = polymec_malloc(c->num_indices*msg->data_size*msg->stride);
    msg->source_procs[i] = proc;
    ++i;
  }
  msg->requests = polymec_malloc((num_sends+num_receives)*sizeof(MPI_Request));
}
#undef HANDLE_PACKING

#define HANDLE_UNPACKING(msg, mpi_type, c_type) \
  if (msg->type == mpi_type) \
  { \
    c_type* src = msg->receive_buffers[i]; \
    c_type* dest = data; \
    for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j) \
      for (int s = 0; s < stride; ++s) \
        dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s]; \
  }

static void mpi_message_unpack(mpi_message_t* msg, 
                               void* data, 
                               ssize_t receive_offset,
                               exchanger_map_t* receive_map)
{
  int pos = 0, proc, i = 0;
  exchanger_channel_t* c;
  while (exchanger_map_next(receive_map, &pos, &proc, &c))
  {
    int* recv_indices = c->indices;
    int stride = msg->stride;
    HANDLE_UNPACKING(msg, MPI_REAL_T, real_t)
    HANDLE_UNPACKING(msg, MPI_FLOAT, float)
    HANDLE_UNPACKING(msg, MPI_INT, int)
    HANDLE_UNPACKING(msg, MPI_LONG, long)
    HANDLE_UNPACKING(msg, MPI_LONG_LONG, long long)
    HANDLE_UNPACKING(msg, MPI_UINT64_T, uint64_t)
    HANDLE_UNPACKING(msg, MPI_INT64_T, int64_t)
    HANDLE_UNPACKING(msg, MPI_CHAR, char)
    HANDLE_UNPACKING(msg, MPI_BYTE, uint8_t)
    ++i;
  }
}
#undef HANDLE_UNPACKING

static void mpi_message_free(mpi_message_t* msg)
{
  if (msg->send_buffers != NULL)
  {
    for (int i = 0; i < msg->num_sends; ++i)
    {
      if (msg->send_buffers[i] != NULL)
        polymec_free(msg->send_buffers[i]);
    }
    polymec_free(msg->send_buffers);
  }
  if (msg->send_buffer_sizes != NULL)
    polymec_free(msg->send_buffer_sizes);
  if (msg->dest_procs != NULL)
    polymec_free(msg->dest_procs);
  if (msg->receive_buffers != NULL)
  {
    for (int i = 0; i < msg->num_receives; ++i)
    {
      if (msg->receive_buffers[i] != NULL)
        polymec_free(msg->receive_buffers[i]);
    }
    polymec_free(msg->receive_buffers);
  }
  if (msg->receive_buffer_sizes != NULL)
    polymec_free(msg->receive_buffer_sizes);
  if (msg->source_procs != NULL)
    polymec_free(msg->source_procs);
  if (msg->requests != NULL)
    polymec_free(msg->requests);
  polymec_free(msg);
}

#if 0
static void mpi_message_fprintf(mpi_message_t* msg, FILE* stream)
{
  if (stream == NULL) return;
  char typeStr[1024];
  if (msg->type == MPI_REAL_T)
    strcpy(typeStr, "real");
  else if (msg->type == MPI_DOUBLE)
    strcpy(typeStr, "double");
  else if (msg->type == MPI_FLOAT)
    strcpy(typeStr, "float");
  else if (msg->type == MPI_INT)
    strcpy(typeStr, "int");
  else if (msg->type == MPI_LONG)
    strcpy(typeStr, "long");
  else if (msg->type == MPI_LONG_LONG)
    strcpy(typeStr, "long long");
  else if (msg->type == MPI_UINT64_T)
    strcpy(typeStr, "uint64_t");
  else
    strcpy(typeStr, "char");
  fprintf(stream, "Message(datatype = %s, stride = %d, tag = %d):\n",
          typeStr, msg->stride, msg->tag);
  if (msg->num_sends > 0)
    fprintf(stream, "Send buffers:\n");
  for (int i = 0; i < msg->num_sends; ++i)
    fprintf(stream, " %d: %d items\n", i, msg->send_buffer_sizes[i]);
  if (msg->num_receives > 0)
    fprintf(stream, "Receive buffers:\n");
  for (int i = 0; i < msg->num_receives; ++i)
    fprintf(stream, " %d: %d items\n", i, msg->receive_buffer_sizes[i]);
}
#endif

#endif

struct exchanger_t
{
  MPI_Comm comm;
  int rank, nprocs;

  // Communication maps.
  exchanger_map_t* send_map;
  exchanger_map_t* receive_map;

  int max_send, max_receive;

  // Send, receive offsets.
  ssize_t send_offset, receive_offset;

  // Pending messages.
  int num_pending_msgs;
  int pending_msg_cap;
  mpi_message_t** pending_msgs;
  void** orig_buffers;
  int** transfer_counts;

  // Deadlock detection.
  real_t dl_thresh;
  int dl_output_rank;
  FILE* dl_output_stream;
};

static void exchanger_clear(exchanger_t* ex)
{
  exchanger_map_clear(ex->send_map);
  exchanger_map_clear(ex->receive_map);

  polymec_free(ex->pending_msgs);
  polymec_free(ex->orig_buffers);
  polymec_free(ex->transfer_counts);

  ex->max_send = -1;
  ex->max_receive = -1;
}

static void exchanger_free(void* ctx)
{
  exchanger_t* ex = ctx;
  exchanger_clear(ex);
  exchanger_map_free(ex->send_map);
  exchanger_map_free(ex->receive_map);
}

exchanger_t* exchanger_new_with_rank(MPI_Comm comm, int rank)
{
  exchanger_t* ex = polymec_gc_malloc(sizeof(exchanger_t), exchanger_free);
  ex->comm = comm;
  ex->rank = rank;
  MPI_Comm_size(comm, &(ex->nprocs));
  ex->send_offset = 0;
  ex->receive_offset = 0;
  ex->dl_thresh = -1.0;
  ex->dl_output_rank = -1;
  ex->dl_output_stream = NULL;
  ex->send_map = exchanger_map_new();
  ex->receive_map = exchanger_map_new();
  ex->num_pending_msgs = 0;
  ex->pending_msg_cap = 32;
  ex->pending_msgs = polymec_malloc(ex->pending_msg_cap * sizeof(mpi_message_t*));
  memset(ex->pending_msgs, 0, ex->pending_msg_cap * sizeof(mpi_message_t*));
  ex->orig_buffers = polymec_malloc(ex->pending_msg_cap * sizeof(void*));
  memset(ex->orig_buffers, 0, ex->pending_msg_cap * sizeof(void*));
  ex->transfer_counts = polymec_malloc(ex->pending_msg_cap * sizeof(int*));
  memset(ex->transfer_counts, 0, ex->pending_msg_cap * sizeof(int*));
  ex->max_send = -1;
  ex->max_receive = -1;

  return ex;
}

exchanger_t* exchanger_new(MPI_Comm comm)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  return exchanger_new_with_rank(comm, rank);
}

exchanger_t* exchanger_clone(exchanger_t* ex)
{
  exchanger_t* clone = exchanger_new(ex->comm);
  int pos = 0, proc;
  int *indices, num_indices;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
    exchanger_set_send(clone, proc, indices, num_indices, true);   
  pos = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
    exchanger_set_receive(clone, proc, indices, num_indices, true);
  return clone;
}

MPI_Comm exchanger_comm(exchanger_t* ex)
{
  return ex->comm;
}

static void delete_map_entry(int key, exchanger_channel_t* value)
{
  exchanger_channel_free(value);
}

void exchanger_set_send(exchanger_t* ex, int remote_process, int* indices, int num_indices, bool copy_indices)
{
  ASSERT(remote_process >= 0);
  ASSERT(remote_process != ex->rank);
  ASSERT(remote_process < ex->nprocs);
  exchanger_channel_t* c = exchanger_channel_new(num_indices, indices, copy_indices);
  exchanger_map_insert_with_kv_dtor(ex->send_map, remote_process, c, delete_map_entry);

  if (remote_process > ex->max_send)
    ex->max_send = remote_process;
}

void exchanger_set_sends(exchanger_t* ex, int_ptr_unordered_map_t* send_map)
{
  int pos = 0, send_proc;
  int_array_t* send_indices;
  while (int_ptr_unordered_map_next(send_map, &pos, &send_proc, (void**)&send_indices))
    exchanger_set_send(ex, send_proc, send_indices->data, (int)send_indices->size, true);   
}

void exchanger_set_send_offset(exchanger_t* ex, ssize_t offset)
{
  ex->send_offset = offset;
}

int exchanger_num_sends(exchanger_t* ex)
{
  return ex->send_map->size;
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

bool exchanger_next_send(exchanger_t* ex, int* pos, int* remote_process, int** indices, int* num_indices)
{
  exchanger_channel_t* c;
  bool result = exchanger_map_next(ex->send_map, pos, remote_process, &c);
  if (result)
  {
    *indices = c->indices;
    *num_indices = c->num_indices;
  }
  return result;
}

bool exchanger_get_send(exchanger_t* ex, int remote_process, int** indices, int* num_indices)
{
  exchanger_channel_t** c_ptr = exchanger_map_get(ex->send_map, remote_process);
  if (c_ptr == NULL)
    return false;
  else
  {
    exchanger_channel_t* c = *c_ptr;
    *indices = c->indices;
    *num_indices = c->num_indices;
    return true;
  }
}

void exchanger_set_receive(exchanger_t* ex, int remote_process, int* indices, int num_indices, bool copy_indices)
{
  ASSERT(remote_process >= 0);
  ASSERT(remote_process != ex->rank);
  ASSERT(remote_process < ex->nprocs);
  exchanger_channel_t* c = exchanger_channel_new(num_indices, indices, copy_indices);
  exchanger_map_insert_with_kv_dtor(ex->receive_map, remote_process, c, delete_map_entry);

  if (remote_process > ex->max_receive)
    ex->max_receive = remote_process;
}

void exchanger_set_receives(exchanger_t* ex, int_ptr_unordered_map_t* recv_map)
{
  int pos = 0, recv_proc;
  int_array_t* recv_indices;
  while (int_ptr_unordered_map_next(recv_map, &pos, &recv_proc, (void**)&recv_indices))
    exchanger_set_receive(ex, recv_proc, recv_indices->data, (int)recv_indices->size, true);   
}

void exchanger_set_receive_offset(exchanger_t* ex, ssize_t offset)
{
  ex->receive_offset = offset;
}

int exchanger_num_receives(exchanger_t* ex)
{
  return ex->receive_map->size;
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

bool exchanger_next_receive(exchanger_t* ex, int* pos, int* remote_process, int** indices, int* num_indices)
{
  exchanger_channel_t* c;
  bool result = exchanger_map_next(ex->receive_map, pos, remote_process, &c);
  if (result)
  {
    *indices = c->indices;
    *num_indices = c->num_indices;
  }
  return result;
}

bool exchanger_get_receive(exchanger_t* ex, int remote_process, int** indices, int* num_indices)
{
  exchanger_channel_t** c_ptr = exchanger_map_get(ex->receive_map, remote_process);
  if (c_ptr == NULL)
    return false;
  else
  {
    exchanger_channel_t* c = *c_ptr;
    *indices = c->indices;
    *num_indices = c->num_indices;
    return true;
  }
}

bool exchanger_verify(exchanger_t* ex, void (*handler)(const char* format, ...))
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();

  // An exchanger is valid/consistent iff the number of elements that 
  // are exchanged between any two processors are agreed upon between those 
  // two processors. 
  
  // First question: do our neighbors agree with us about being our 
  // neighbors?

  // Tally up the number of indices we're sending to and receiving from 
  // every other process in the communicator.
  int pos = 0, proc, *indices, num_indices;
  int my_neighbors[2*ex->nprocs];
  memset(my_neighbors, 0, 2 * sizeof(int) * ex->nprocs);
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
    my_neighbors[2*proc] += num_indices;
  pos = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
    my_neighbors[2*proc+1] += num_indices;

  // Do an all-to-all exchange to get everyone's votes on who is whose 
  // neighbor.
  int* neighbors_for_proc = polymec_malloc(sizeof(int)*2*ex->nprocs*ex->nprocs);
  MPI_Allgather(my_neighbors, 2*ex->nprocs, MPI_INT, 
                neighbors_for_proc, 2*ex->nprocs, MPI_INT, ex->comm);

set_log_mode(LOG_TO_ALL_RANKS);
set_log_level(LOG_DEBUG);
  for (int p = 0; p < ex->nprocs; ++p)
  {
    if (p != ex->rank)
    {
      int num_im_sending = my_neighbors[2*p];
      int num_theyre_receiving = neighbors_for_proc[2*(p*ex->nprocs+ex->rank)+1];
      if (num_im_sending != num_theyre_receiving)
      {
        polymec_free(neighbors_for_proc);
        if (handler != NULL)
          handler("exchanger_verify: Proc %d is sending %d elements to proc %d, which is expecting %d elements.",
              ex->rank, num_im_sending, p, num_theyre_receiving);
        return false;
      }

      int num_im_receiving = my_neighbors[2*p+1];
      int num_theyre_sending = neighbors_for_proc[2*(p*ex->nprocs+ex->rank)];
      if (num_im_receiving != num_theyre_sending)
      {
        polymec_free(neighbors_for_proc);
        if (handler != NULL)
          handler("exchanger_verify: Proc %d is sending %d elements to proc %d, which is expecting %d elements.",
              p, num_theyre_sending, ex->rank, num_im_receiving);
        return false;
      }
    }
  }
  polymec_free(neighbors_for_proc);

  // Next: do our indices agree with those of our neighbors?

  // So we send our expected number of exchanged elements 
  // to our neighbors, and receive their numbers, and then compare.
  int num_send_procs = ex->send_map->size;
  int num_receive_procs = ex->receive_map->size;
  int num_sent_elements[num_send_procs], 
      num_elements_expected_by_neighbors[num_send_procs],
      send_procs[num_send_procs];
  MPI_Request requests[num_send_procs + num_receive_procs];

  // Post receives for numbers of elements we send to others.
  int p = 0;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    num_sent_elements[p] = num_indices;
    send_procs[p] = proc;
    int err = MPI_Irecv(&num_elements_expected_by_neighbors[p], 1, 
                        MPI_INT, proc, 0, ex->comm, &requests[p]);
    if (err != MPI_SUCCESS)
    {
      if (handler != NULL)
        handler("exchanger_verify: Could not receive sent element data from %d", proc);
      return false;
    }
    ++p;
  }

  // Send our receive neighbors the number that we expect to receive. This 
  // means that we will get the number that our neighbors expect to receive,
  // which should be the same number as what we're sending.
  pos = p = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
  {
    int err = MPI_Isend(&num_indices, 1, MPI_INT, proc, 0, ex->comm, &requests[num_send_procs + p]);
    if (err != MPI_SUCCESS)
    {
      if (handler != NULL)
        handler("exchanger_verify: Could not send number of received elements to %d", proc);
      return false;
    }
    ++p;
  }

  // Wait for messages to fly.
  MPI_Status statuses[num_send_procs + num_receive_procs];
  MPI_Waitall(num_send_procs + num_receive_procs, requests, statuses);

  // Now verify that the numbers are the same.
  for (p = 0; p < num_send_procs; ++p)
  {
    if (num_sent_elements[p] != num_elements_expected_by_neighbors[p])
    {
      if (handler != NULL)
      {
        handler("exchanger_verify: Sending %d elements to proc %d, but %d elements are expected.", 
                num_sent_elements[p], send_procs[p],
                num_elements_expected_by_neighbors[p]);
      }
      return false;
    }
  }
  STOP_FUNCTION_TIMER();
#endif
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
  START_FUNCTION_TIMER();
  int token = exchanger_start_exchange(ex, data, stride, tag, type);
  exchanger_finish_exchange(ex, token);
  STOP_FUNCTION_TIMER();
}

#if POLYMEC_HAVE_MPI
static int exchanger_send_message(exchanger_t* ex, void* data, mpi_message_t* msg)
{
  START_FUNCTION_TIMER();
  // If we are expecting data, post asynchronous receives. 
  int j = 0;
  for (int i = 0; i < msg->num_receives; ++i)
  {
    ASSERT(ex->rank != msg->source_procs[i]);
    int err = MPI_Irecv(msg->receive_buffers[i], 
                        msg->stride * msg->receive_buffer_sizes[i],
                        msg->type, msg->source_procs[i], msg->tag, ex->comm, 
                        &(msg->requests[j++]));
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error posting receive from %d: %d\n(%s)\n", 
          ex->rank, msg->source_procs[i], err, str);
      polymec_error(err_msg);
    }
  }

  // Send the data asynchronously. 
  for (int i = 0; i < msg->num_sends; ++i)
  {
    int err = MPI_Isend(msg->send_buffers[i], 
                        msg->stride * msg->send_buffer_sizes[i], 
                        msg->type, msg->dest_procs[i], msg->tag, ex->comm, 
                        &(msg->requests[j++])); 
    if (err != MPI_SUCCESS)
    {
      int resultlen;
      char str[MPI_MAX_ERROR_STRING];
      MPI_Error_string(err, str, &resultlen);
      char err_msg[1024];
      snprintf(err_msg, 1024, "%d: MPI Error sending to %d: %d\n(%s)\n", 
          ex->rank, msg->dest_procs[i], err, str);
      polymec_error(err_msg);
    }
  }
  msg->num_requests = j;

  // Allocate a token.
  int token = 0;
  while ((token < ex->num_pending_msgs) && (ex->pending_msgs[token] != NULL))
    ++token;
  if (token == ex->num_pending_msgs) ++ex->num_pending_msgs;
  if (ex->num_pending_msgs == ex->pending_msg_cap)
  {
    ex->pending_msg_cap *= 2;
    ex->pending_msgs = polymec_realloc(ex->pending_msgs, ex->pending_msg_cap*sizeof(mpi_message_t));
    ex->orig_buffers = polymec_realloc(ex->orig_buffers, ex->pending_msg_cap*sizeof(void*));
    ex->transfer_counts = polymec_realloc(ex->transfer_counts, ex->pending_msg_cap*sizeof(int*));
  }
  ex->pending_msgs[token] = msg;
  ex->orig_buffers[token] = data;
  STOP_FUNCTION_TIMER();
  return token;
}
#endif

int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  // Create a message for this array.
  mpi_message_t* msg = mpi_message_new(type, stride, tag);
  mpi_message_pack(msg, data, ex->send_offset, ex->send_map, ex->receive_map);
  
  // Begin the transmission and allocate a token for it.
  int token = exchanger_send_message(ex, data, msg);
  STOP_FUNCTION_TIMER();
  return token;
#else
  return 0;
#endif
}

#if POLYMEC_HAVE_MPI
static int exchanger_waitall(exchanger_t* ex, mpi_message_t* msg)
{
  // Allocate storage for statuses of sends/receives.
  int num_requests = msg->num_requests;
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
          if (MPI_Test(&(msg->requests[i]), &(finished[i]), &(statuses[i])) != MPI_SUCCESS)
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
            MPI_Cancel(&(msg->requests[i]));
        }

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
            if (expecting_data && (i < ex->receive_map->size)) // outstanding receive 
            {
              outstanding_receive_procs[num_outstanding_receives] = msg->source_procs[i];
              outstanding_receive_bytes[num_outstanding_receives] = msg->receive_buffer_sizes[i] * msg->data_size;
              ++num_outstanding_receives;
            }
            else if (sent_data) // outstanding send 
            {
              outstanding_send_procs[num_outstanding_sends] = msg->dest_procs[i - msg->num_receives];
              outstanding_send_bytes[num_outstanding_sends] = msg->send_buffer_sizes[i - msg->num_receives] * msg->data_size;
              ++num_outstanding_sends;
            }
          }
          else
          {
            if (expecting_data && (i < msg->num_receives)) // completed receive 
            {
              completed_receive_procs[num_completed_receives] = msg->source_procs[i];
              completed_receive_bytes[num_completed_receives] = msg->receive_buffer_sizes[i] * msg->data_size;
              ++num_completed_receives;
            }
            else if (sent_data) // completed send 
            {
              completed_send_procs[num_completed_sends] = msg->dest_procs[i - msg->num_receives];
              completed_send_bytes[num_completed_sends] = msg->send_buffer_sizes[i - msg->num_receives] * msg->data_size;
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
        if (i < msg->num_receives)
        {
          // Now we can really get nitty-gritty and try to diagnose the
          // problem carefully! 
          if (statuses[i].MPI_ERROR == MPI_ERR_TRUNCATE)
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n"
                    "(Expected %d bytes)\n", ex->rank, msg->source_procs[i], statuses[i].MPI_ERROR, 
                    errstr, msg->receive_buffer_sizes[i]);
          }
          else
          {
            fprintf(ex->dl_output_stream, "%d: MPI error receiving from %d (%d) %s\n",
                    ex->rank, msg->source_procs[i], statuses[i].MPI_ERROR, errstr);
          }
        }
        else 
        {
          fprintf(ex->dl_output_stream, "%d: MPI error sending to %d (%d) %s\n",
                  ex->rank, msg->dest_procs[i - msg->num_receives], statuses[i].MPI_ERROR, errstr);
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
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);
  ASSERT(ex->transfer_counts[token] == NULL); // We can't finish a transfer as an exchange!

  // Retrieve the message for the given token.
  mpi_message_t* msg = ex->pending_msgs[token];
  int err = exchanger_waitall(ex, msg);
  if (err != MPI_SUCCESS) return;

  // Copy the received data into the original array.
  void* orig_buffer = ex->orig_buffers[token];
  mpi_message_unpack(msg, orig_buffer, ex->receive_offset, ex->receive_map);

  // Pull the message out of our list of pending messages and delete it.
  ex->pending_msgs[token] = NULL;
  ex->orig_buffers[token] = NULL;
  mpi_message_free(msg);
  STOP_FUNCTION_TIMER();
#endif
}

void
exchanger_fprintf(exchanger_t* ex, FILE* stream)
{
#if POLYMEC_HAVE_MPI
  fprintf(stream, "Exchanger(");
  if (ex->comm == MPI_COMM_WORLD)
    fprintf(stream, "MPI_COMM_WORLD");
  else
    fprintf(stream, "[user communicator]");
  fprintf(stream, ", rank %d):", ex->rank);
  if ((ex->send_map->size == 0) && (ex->receive_map->size == 0))
    fprintf(stream, " (no transactions)\n");
  else
  {
    fprintf(stream, "\n");
    int pos = 0, proc;
    exchanger_channel_t* c;
    while (exchanger_map_next(ex->send_map, &pos, &proc, &c))
    {
      fprintf(stream, " [%d] -> [%d]: ", ex->rank, proc);
      for (int j = 0; j < c->num_indices; ++j)
        fprintf(stream, " %d ", c->indices[j]);
      fprintf(stream, " (%d elements)\n", c->num_indices);
    }
    pos = 0;
    while (exchanger_map_next(ex->receive_map, &pos, &proc, &c))
    {
      fprintf(stream, " [%d] <- [%d]: ", ex->rank, proc);
      for (int j = 0; j < c->num_indices; ++j)
        fprintf(stream, " %d ", c->indices[j]);
      fprintf(stream, " (%d elements)\n", c->num_indices);
    }
  }
#else
  fprintf(stream, "Exchanger(dummy)\n");
#endif
}

static size_t ex_size(void* obj)
{
#if POLYMEC_HAVE_MPI
  exchanger_t* ex = obj;
  size_t size = 2 * sizeof(int);
  int pos = 0, proc, *indices, num_indices;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
    size += (2 + num_indices) * sizeof(int);
  pos = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
    size += (2 + num_indices) * sizeof(int);
  return size;
#else
  return 0;
#endif
}

static void* ex_read(byte_array_t* bytes, size_t* offset)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  exchanger_t* ex = exchanger_new(comm);
#if POLYMEC_HAVE_MPI
  int num_sends, num_receives;
  MPI_Comm_rank(comm, &ex->rank);
  byte_array_read_ints(bytes, 1, &num_sends, offset);
  for (int i = 0; i < num_sends; ++i)
  {
    int proc, num_indices;
    byte_array_read_ints(bytes, 1, &proc, offset);
    byte_array_read_ints(bytes, 1, &num_indices, offset);
    int indices[num_indices];
    byte_array_read_ints(bytes, num_indices, indices, offset);
    exchanger_set_send(ex, proc, indices, num_indices, true);
  }
  byte_array_read_ints(bytes, 1, &num_receives, offset);
  for (int i = 0; i < num_receives; ++i)
  {
    int proc, num_indices;
    byte_array_read_ints(bytes, 1, &proc, offset);
    byte_array_read_ints(bytes, 1, &num_indices, offset);
    int indices[num_indices];
    byte_array_read_ints(bytes, num_indices, indices, offset);
    exchanger_set_receive(ex, proc, indices, num_indices, true);
  }
#endif 
  return ex;
}

static void ex_write(void* obj, byte_array_t* bytes, size_t* offset)
{
#if POLYMEC_HAVE_MPI
  exchanger_t* ex = obj;

  int pos = 0, proc, *indices, num_indices;
  byte_array_write_ints(bytes, 1, &ex->send_map->size, offset);
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    byte_array_write_ints(bytes, 1, &proc, offset);
    byte_array_write_ints(bytes, 1, &num_indices, offset);
    byte_array_write_ints(bytes, num_indices, indices, offset);
  }
  byte_array_write_ints(bytes, 1, &ex->receive_map->size, offset);
  pos = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
  {
    byte_array_write_ints(bytes, 1, &proc, offset);
    byte_array_write_ints(bytes, 1, &num_indices, offset);
    byte_array_write_ints(bytes, num_indices, indices, offset);
  }
#endif
}

serializer_t* exchanger_serializer()
{
  return serializer_new("exchanger", ex_size, ex_read, ex_write, NULL);
}

void** exchanger_create_metadata_send_arrays(exchanger_t* ex,
                                             MPI_Datatype type,
                                             int stride)
{
#if POLYMEC_HAVE_MPI
  int num_arrays = exchanger_num_sends(ex);
  uint8_t** arrays = polymec_malloc(sizeof(void*) * num_arrays);
  int pos = 0, i = 0, proc, *indices, num_indices;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
    arrays[i] = polymec_malloc(mpi_size(type) * num_indices);
  return (void**)arrays;
#else
  return NULL;
#endif
}

void exchanger_free_metadata_send_arrays(exchanger_t* ex, 
                                         void** arrays)
{
#if POLYMEC_HAVE_MPI
  int num_arrays = exchanger_num_sends(ex);
  for (int i = 0; i < num_arrays; ++i)
    polymec_free(arrays[i]);
  polymec_free(arrays);
#endif
}

void** exchanger_create_metadata_receive_arrays(exchanger_t* ex,
                                                MPI_Datatype type,
                                                int stride)
{
#if POLYMEC_HAVE_MPI
  int num_arrays = exchanger_num_receives(ex);
  uint8_t** arrays = polymec_malloc(sizeof(void*) * num_arrays);
  int pos = 0, i = 0, proc, *indices, num_indices;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
    arrays[i] = polymec_malloc(mpi_size(type) * num_indices);
  return (void**)arrays;
#else
  return NULL;
#endif
}

void exchanger_free_metadata_receive_arrays(exchanger_t* ex, 
                                            void** arrays)
{
#if POLYMEC_HAVE_MPI
  int num_arrays = exchanger_num_receives(ex);
  for (int i = 0; i < num_arrays; ++i)
    polymec_free(arrays[i]);
  polymec_free(arrays);
#endif
}

int exchanger_start_metadata_transfer(exchanger_t* ex,
                                      void** send_arrays,
                                      void** receive_arrays,
                                      int stride,
                                      int tag,
                                      MPI_Datatype type,
                                      exchanger_metadata_dir direction)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();

  // Create a message using these metadata arrays, and set it up manually.
  mpi_message_t* msg = mpi_message_new(type, stride, tag);

  int num_sends = exchanger_num_sends(ex);
  int* dest_procs = polymec_malloc(sizeof(int) * num_sends);
  int* send_buffer_sizes = polymec_malloc(sizeof(int) * num_sends);
  int pos = 0, proc, *indices, num_indices, i = 0;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    dest_procs[i] = proc;
    send_buffer_sizes[i++] = num_indices;
  }

  int num_receives = exchanger_num_receives(ex);
  int* source_procs = polymec_malloc(sizeof(int) * num_receives);
  int* receive_buffer_sizes = polymec_malloc(sizeof(int) * num_receives);
  pos = 0; i = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
  {
    source_procs[i] = proc;
    receive_buffer_sizes[i++] = num_indices;
  }

  // Set things up based on the direction.
  if (direction == EX_METADATA_FORWARD)
  {
    msg->num_sends = num_sends;
    msg->dest_procs = dest_procs;
    msg->send_buffer_sizes = send_buffer_sizes;
    msg->send_buffers = send_arrays;
    msg->source_procs = source_procs;
    msg->receive_buffer_sizes = receive_buffer_sizes;
    msg->receive_buffers = receive_arrays;
  }
  else
  {
    // Reverse!
    msg->num_sends = num_receives;
    msg->dest_procs = source_procs;
    msg->send_buffer_sizes = receive_buffer_sizes;
    msg->send_buffers = receive_arrays;
    msg->source_procs = dest_procs;
    msg->receive_buffer_sizes = send_buffer_sizes;
    msg->receive_buffers = send_arrays;
  }

  // Allocate requests.
  msg->requests = polymec_malloc((msg->num_sends+msg->num_receives)*sizeof(MPI_Request));

  int token = exchanger_send_message(ex, NULL, msg);
  STOP_FUNCTION_TIMER();
  return token;
#else
  return 0;
#endif
}

void exchanger_finish_metadata_transfer(exchanger_t* ex,
                                        int token)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);
  ASSERT(ex->transfer_counts[token] == NULL); // We can't finish a transfer as an exchange!
  ASSERT(ex->orig_buffers[token] == NULL); // This is a metadata transfer, right?

  // Retrieve the message for the given token.
  mpi_message_t* msg = ex->pending_msgs[token];
  int err = exchanger_waitall(ex, msg);
  if (err != MPI_SUCCESS) return;

  // Release the send/receive buffers.
  msg->send_buffers = NULL;
  msg->receive_buffers = NULL;

  // Pull the message out of our list of pending messages and delete it.
  ex->pending_msgs[token] = NULL;
  mpi_message_free(msg);
  STOP_FUNCTION_TIMER();
#endif
}

void exchanger_transfer_metadata(exchanger_t* ex,
                                 void** send_arrays,
                                 void** receive_arrays,
                                 int stride,
                                 int tag,
                                 MPI_Datatype type,
                                 exchanger_metadata_dir direction)
{
  int token = exchanger_start_metadata_transfer(ex, send_arrays, receive_arrays, stride, tag, type, direction);
  exchanger_finish_metadata_transfer(ex, token);
}

struct migrator_t
{
  exchanger_t* ex;
};

static void migrator_free(void* ctx)
{
  migrator_t* m = ctx;
  m->ex = NULL;
}

migrator_t* migrator_new(MPI_Comm comm)
{
  migrator_t* m = polymec_gc_malloc(sizeof(migrator_t), migrator_free);
  m->ex = exchanger_new(comm);
  return m;
}

migrator_t* migrator_clone(migrator_t* m)
{
  migrator_t* clone = polymec_gc_malloc(sizeof(migrator_t), migrator_free);
  clone->ex = exchanger_clone(m->ex);
  return clone;
}

int migrator_num_sends(migrator_t* m)
{
  return exchanger_num_sends(m->ex);
}

bool migrator_next_send(migrator_t* m, int* pos, int* remote_process, int** indices, int* num_indices)
{
  return exchanger_next_send(m->ex, pos, remote_process, indices, num_indices);
}

int migrator_num_receives(migrator_t* m)
{
  return exchanger_num_receives(m->ex);
}

bool migrator_next_receive(migrator_t* m, int* pos, int* remote_process, int** indices, int* num_indices)
{
  return exchanger_next_receive(m->ex, pos, remote_process, indices, num_indices);
}

void migrator_verify(migrator_t* m, void (*handler)(const char* format, ...))
{
  exchanger_verify(m->ex, handler);
}

migrator_t* migrator_from_global_partition(MPI_Comm comm, 
                                           int64_t* global_partition,
                                           int num_global_vertices)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  exchanger_t* ex = exchanger_new(comm);
  
  if (rank == 0)
  {
    // Get the number of vertices we're going to send to each other process.
    int num_vertices_to_send[nprocs];
    memset(num_vertices_to_send, 0, sizeof(int) * nprocs);
    for (int v = 0; v < num_global_vertices; ++v)
      num_vertices_to_send[global_partition[v]]++;

    // Send this number.
    int n;
    MPI_Scatter(num_vertices_to_send, 1, MPI_INT, &n, 1, MPI_INT, 0, comm);

    // Now send the vertices to each process and register these sends
    // with the distributor.
    for (int p = 1; p < nprocs; ++p)
    {
      ASSERT(num_vertices_to_send[p] > 0);
      int vertices[num_vertices_to_send[p]], k = 0;
      for (int i = 0; i < num_global_vertices; ++i)
      {
        if (global_partition[i] == p)
          vertices[k++] = i;
      }
      exchanger_set_send(ex, p, vertices, num_vertices_to_send[p], true);
    }

    // Figure out the local vertices for rank 0.
    int num_local_vertices = num_vertices_to_send[0];
    int local_vertices[num_local_vertices], k = 0;
    for (int i = 0; i < num_global_vertices; ++i)
    {
      if (global_partition[i] == 0)
        local_vertices[k++] = i;
    }
  }
  else
  {
    // Get the number of vertices we will receive from rank 0.
    int num_local_vertices;
    MPI_Scatter(NULL, 1, MPI_INT, &num_local_vertices, 1, MPI_INT, 0, comm);

    // Now get the vertices.
    int local_vertices[num_local_vertices];
    for (int i = 0; i < num_local_vertices; ++i)
      local_vertices[i] = i;

    // Now register all the vertices we're receiving with the distributor.
    ASSERT(num_local_vertices > 0);
    exchanger_set_receive(ex, 0, local_vertices, num_local_vertices, true);
  }

  migrator_t* m = polymec_gc_malloc(sizeof(migrator_t), migrator_free);
  m->ex = ex;
  STOP_FUNCTION_TIMER();
  return m;
}

migrator_t* migrator_from_local_partition(MPI_Comm comm,
                                          int64_t* local_partition,
                                          int num_vertices)
{
  START_FUNCTION_TIMER();
  exchanger_t* ex = exchanger_new(comm);
  
  // Tally up the number of vertices we're going to send to every other process, 
  // including those that are staying on our own.
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  int num_vertices_to_send[nprocs];
  memset(num_vertices_to_send, 0, sizeof(int) * nprocs);
  for (int v = 0; v < num_vertices; ++v)
    num_vertices_to_send[local_partition[v]]++;

  // Get the number of vertices we're going to receive from every other process.
  int num_vertices_to_receive[nprocs];
  MPI_Alltoall(num_vertices_to_send, 1, MPI_INT, 
               num_vertices_to_receive, 1, MPI_INT, comm);

  // Set up send and receive vertices.
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
          int receive_vertices[num_vertices_to_receive[p]];
          for (int i = 0; i < num_vert; ++i, ++recv_offset)
            receive_vertices[i] = num_local_vert + recv_offset;
          exchanger_set_receive(ex, p, receive_vertices, num_vertices_to_receive[p], true);
        }
      }
      if (num_vertices_to_send[p] > 0)
      {
        int send_vertices[num_vertices_to_send[p]], s = 0;
        for (int v = 0; v < num_vertices; ++v)
        {
          if (local_partition[v] == p)
            send_vertices[s++] = v;
        }
        if (p != rank)
          exchanger_set_send(ex, p, send_vertices, num_vertices_to_send[p], true);
      }
    }
  }

  migrator_t* m = polymec_gc_malloc(sizeof(migrator_t), migrator_free);
  m->ex = ex;
  STOP_FUNCTION_TIMER();
  return m;
}

void migrator_transfer(migrator_t* m, void* data, int* count, int stride, int tag, MPI_Datatype type)
{
  START_FUNCTION_TIMER();
  int token = migrator_start_transfer(m, data, count, stride, tag, type);
  migrator_finish_transfer(m, token);
  STOP_FUNCTION_TIMER();
}

int migrator_start_transfer(migrator_t* m, void* data, int* count, int stride, int tag, MPI_Datatype type)
{
#if POLYMEC_HAVE_MPI
  // This is the same as an exchange, with an additional datum (the length 
  // of data) stored.
  int token = exchanger_start_exchange(m->ex, data, stride, tag, type);
  m->ex->transfer_counts[token] = count; 
  return token;
#else
  return 0;
#endif
}

void migrator_finish_transfer(migrator_t* m, int token)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT(token < m->ex->num_pending_msgs);
  ASSERT(m->ex->transfer_counts[token] != NULL); // We can't finish an exchange as a transfer!

  // Retrieve the message and its size (in data elements) for the given token.
  mpi_message_t* msg = m->ex->pending_msgs[token];
  int* count = m->ex->transfer_counts[token];

  // Wait for the transfer to finish.
  int err = exchanger_waitall(m->ex, msg);
  if (err != MPI_SUCCESS) return;

  // Copy the received data into the original array.
  void* orig_buffer = m->ex->orig_buffers[token];
  mpi_message_unpack(msg, orig_buffer, m->ex->receive_offset, m->ex->receive_map);

  // Count up the elements we actually received.
  int num_received = 0, pos = 0, proc;
  exchanger_channel_t* c;
  while (exchanger_map_next(m->ex->receive_map, &pos, &proc, &c))
    num_received += c->num_indices;

  // We need to kill all of the sent elements that aren't overwritten by received elements.
  // So first make a list of all the received indices.
  int num_sent = 0;
  int_unordered_set_t* overwritten_sends = int_unordered_set_new();
  pos = 0;
  while (exchanger_map_next(m->ex->send_map, &pos, &proc, &c))
  {
    for (int i = 0; i < c->num_indices; ++i)
      int_unordered_set_insert(overwritten_sends, c->indices[i]);
    num_sent += c->num_indices;
  }

  // Update the count.
  *count += (num_received - num_sent);
  ASSERT(*count >= 0);
  if (num_sent - num_received > 0)
  {
    log_debug("migrator_transfer: deleting %d elements locally.", 
              num_sent - num_received);
  }
  else if (num_sent - num_received < 0)
  {
    log_debug("migrator_transfer: adding %d elements locally.", 
              num_received - num_sent);
  }

  // Record the "width" in bytes of the type we're transferring.
  size_t width = mpi_size(msg->type);

  // Remove the send elements that aren't overwritten by received ones.
  // We recreate the array in buffer, skipping over the elements we're removing.
  char* array = orig_buffer;
  int stride = msg->stride, offset = 0;
  char* buffer = polymec_malloc(width * stride * (*count));
  pos = 0;
  while (exchanger_map_next(m->ex->send_map, &pos, &proc, &c))
  {
    for (int i = 0; i < c->num_indices; ++i)
    {
      int index = c->indices[i];
      if (!int_unordered_set_contains(overwritten_sends, index))
      {
        memcpy(&buffer[stride*width*offset], &array[stride*width*index], stride*width);
        ++offset;
      }
    }
  }
  memcpy(array, buffer, width*stride*offset);
  polymec_free(buffer);
  int_unordered_set_free(overwritten_sends);

  // Pull the message out of our list of pending messages and delete it.
  m->ex->pending_msgs[token] = NULL;
  m->ex->orig_buffers[token] = NULL;
  m->ex->transfer_counts[token] = NULL;
  mpi_message_free(msg);
  STOP_FUNCTION_TIMER();
#endif
}

void migrator_fprintf(migrator_t* m, FILE* stream)
{
  fprintf(stream, "Migrator: ");
  exchanger_fprintf(m->ex, stream);
}

static size_t m_size(void* obj)
{
  return ex_size(obj);
}

static void* m_read(byte_array_t* bytes, size_t* offset)
{
  migrator_t* m = polymec_gc_malloc(sizeof(migrator_t), migrator_free);
  m->ex = ex_read(bytes, offset);
  return m;
}

static void m_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  migrator_t* m = obj;
  ex_write(m->ex, bytes, offset);
}

serializer_t* migrator_serializer()
{
  return serializer_new("migrator", m_size, m_read, m_write, NULL);
}

