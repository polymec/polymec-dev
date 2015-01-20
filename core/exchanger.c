// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "exchanger.h"
#include "core/unordered_map.h"
#include "core/array.h"

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
static mpi_message_t* mpi_message_new(MPI_Datatype type, int stride, int tag)
{
  ASSERT(stride > 0);
  mpi_message_t* msg = polymec_malloc(sizeof(mpi_message_t));
  msg->type = type;
  msg->stride = stride;
  msg->tag = tag;
  if (type == MPI_REAL_T)
    msg->data_size = sizeof(real_t);
  else if (type == MPI_DOUBLE)
    msg->data_size = sizeof(double);
  else if (type == MPI_FLOAT)
    msg->data_size = sizeof(float);
  else if (type == MPI_INT)
    msg->data_size = sizeof(int);
  else if (type == MPI_LONG)
    msg->data_size = sizeof(long);
  else if (type == MPI_LONG_LONG)
    msg->data_size = sizeof(long long);
  else if (type == MPI_UINT64_T)
    msg->data_size = sizeof(uint64_t);
  else if (type == MPI_CHAR)
    msg->data_size = sizeof(char);
  else if (type == MPI_BYTE)
    msg->data_size = sizeof(uint8_t);
  else 
    polymec_error("Unsupported MPI data type used to construct MPI message.");
  msg->requests = NULL;
  return msg;
}

static void mpi_message_pack(mpi_message_t* msg, 
                             void* data, 
                             size_t send_offset,
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

    if (msg->type == MPI_REAL_T)
    {
      real_t* src = data;
      real_t* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_DOUBLE)
    {
      double* src = data;
      double* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_FLOAT)
    {
      float* src = data;
      float* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_INT)
    {
      int* src = data;
      int* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_LONG)
    {
      long* src = data;
      long* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_LONG_LONG)
    {
      long long* src = data;
      long long* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_UINT64_T)
    {
      uint64_t* src = data;
      uint64_t* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_CHAR)
    {
      char* src = data;
      char* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    else if (msg->type == MPI_BYTE)
    {
      uint8_t* src = data;
      uint8_t* dest = msg->send_buffers[i];
      for (int j = 0; j < msg->send_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[stride*j+s] = src[send_offset+stride*send_indices[j]+s];
    }
    ++i;
  }

  pos = 0, i = 0;
  while (exchanger_map_next(receive_map, &pos, &proc, &c))
  {
    msg->receive_buffer_sizes[i] = c->num_indices;
    msg->receive_buffers[i] = polymec_malloc(c->num_indices*msg->data_size*msg->stride);
    msg->source_procs[i] = proc;
    ++i;
  }
  msg->requests = polymec_malloc((num_sends+num_receives)*sizeof(MPI_Request));
}

static void mpi_message_unpack(mpi_message_t* msg, 
                               void* data, 
                               size_t receive_offset,
                               exchanger_map_t* receive_map)
{
  int pos = 0, proc, i = 0;
  exchanger_channel_t* c;
  while (exchanger_map_next(receive_map, &pos, &proc, &c))
  {
    int* recv_indices = c->indices;
    int stride = msg->stride;
    if (msg->type == MPI_REAL_T)
    {
      real_t* src = msg->receive_buffers[i];
      real_t* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_DOUBLE)
    {
      double* src = msg->receive_buffers[i];
      double* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_FLOAT)
    {
      float* src = msg->receive_buffers[i];
      float* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_INT)
    {
      int* src = msg->receive_buffers[i];
      int* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_LONG)
    {
      long* src = msg->receive_buffers[i];
      long* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_LONG_LONG)
    {
      long long* src = msg->receive_buffers[i];
      long long* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_UINT64_T)
    {
      uint64_t* src = msg->receive_buffers[i];
      uint64_t* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_CHAR)
    {
      char* src = msg->receive_buffers[i];
      char* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    else if (msg->type == MPI_BYTE)
    {
      uint8_t* src = msg->receive_buffers[i];
      uint8_t* dest = data;
      for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j)
        for (int s = 0; s < stride; ++s)
          dest[receive_offset+stride*recv_indices[j]+s] = src[stride*j+s];
    }
    ++i;
  }
}

static void mpi_message_free(mpi_message_t* msg)
{
  for (int i = 0; i < msg->num_sends; ++i)
  {
    if (msg->send_buffers[i] != NULL)
      polymec_free(msg->send_buffers[i]);
  }
  polymec_free(msg->send_buffers);
  polymec_free(msg->send_buffer_sizes);
  polymec_free(msg->dest_procs);
  for (int i = 0; i < msg->num_receives; ++i)
  {
    if (msg->receive_buffers[i] != NULL)
      polymec_free(msg->receive_buffers[i]);
  }
  polymec_free(msg->receive_buffers);
  polymec_free(msg->receive_buffer_sizes);
  polymec_free(msg->source_procs);
  if (msg->requests != NULL)
    polymec_free(msg->requests);
  polymec_free(msg);
}

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

struct exchanger_t
{
  MPI_Comm comm;
  int rank, nprocs;

  // Communication maps.
  exchanger_map_t* send_map;
  exchanger_map_t* receive_map;

  int max_send, max_receive;

  // Send, receive offsets.
  size_t send_offset, receive_offset;

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

exchanger_t* exchanger_new_with_rank(MPI_Comm comm, int rank)
{
  exchanger_t* ex = polymec_malloc(sizeof(exchanger_t));
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

void exchanger_free(exchanger_t* ex)
{
  exchanger_clear(ex);
  exchanger_map_free(ex->send_map);
  exchanger_map_free(ex->receive_map);
  polymec_free(ex);
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
    exchanger_set_send(ex, send_proc, send_indices->data, send_indices->size, true);   
}

void exchanger_set_send_offset(exchanger_t* ex, size_t offset)
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
    exchanger_set_receive(ex, recv_proc, recv_indices->data, recv_indices->size, true);   
}

void exchanger_set_receive_offset(exchanger_t* ex, size_t offset)
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


void exchanger_verify(exchanger_t* ex, void (*handler)(const char* format, ...))
{
#if POLYMEC_HAVE_MPI
  // An exchanger is valid/consistent iff the number of elements that 
  // are exchanged between any two processors are agreed upon between those 
  // two processors. So we send our expected number of exchanged elements 
  // to our neighbors, and receive their numbers, and then compare.
  int num_send_procs = ex->send_map->size;
  int num_receive_procs = ex->receive_map->size;
  int num_sent_elements[num_send_procs], 
      num_elements_expected_by_neighbors[num_send_procs],
      send_procs[num_send_procs];
  MPI_Request requests[num_send_procs + num_receive_procs];
  int pos = 0, proc, *indices, num_indices;

  // Post receives for numbers of elements we send to others.
  int p = 0;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    num_sent_elements[p] = num_indices;
    send_procs[p] = proc;
    int err = MPI_Irecv(&num_elements_expected_by_neighbors[p], 1, 
                        MPI_INT, proc, 0, ex->comm, &requests[p]);
    if (err != MPI_SUCCESS)
      polymec_error("exchanger_verify: Could not receive sent element data from %d", proc);
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
      polymec_error("exchanger_verify: Could not send number of received elements to %d", proc);
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
      handler("exchanger_verify: Sending %d elements to proc %d, but %d elements are expected.", 
              num_sent_elements[p], send_procs[p],
              num_elements_expected_by_neighbors[p]);
    }
  }
#endif
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
  int token = exchanger_start_exchange(ex, data, stride, tag, type);
  exchanger_finish_exchange(ex, token);
}

int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type)
{
#if POLYMEC_HAVE_MPI
  // Create a message for this array.
  mpi_message_t* msg = mpi_message_new(type, stride, tag);
  mpi_message_pack(msg, data, ex->send_offset, ex->send_map, ex->receive_map);

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

  // Allocate a token for the transmission and store the pending message.
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
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);
  ASSERT(ex->transfer_counts[token] == NULL); // We can't finishing a transfer as an exchange!

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
#endif
}

void exchanger_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type)
{
#if POLYMEC_HAVE_MPI
  int token = exchanger_start_transfer(ex, data, count, stride, tag, type);
  exchanger_finish_transfer(ex, token);
#endif
}

int exchanger_start_transfer(exchanger_t* ex, void* data, int* count, int stride, int tag, MPI_Datatype type)
{
#if POLYMEC_HAVE_MPI
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
#if POLYMEC_HAVE_MPI
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);
  ASSERT(ex->transfer_counts[token] != NULL); // We can't finish an exchange as a transfer!

  // This is like exchanger_finish_exchange, but with some rearranging.

  // Retrieve the message for the given token.
  mpi_message_t* msg = ex->pending_msgs[token];
  int* count = ex->transfer_counts[token];
  int err = exchanger_waitall(ex, msg);
  if (err != MPI_SUCCESS) return;

  // Copy the received data into the original array.
  void* orig_buffer = ex->orig_buffers[token];
  mpi_message_unpack(msg, orig_buffer, ex->receive_offset, ex->receive_map);

  // Now cull the sent elements from the array.
  int num_sent = 0, pos = 0, proc;
  exchanger_channel_t* c;
  while (exchanger_map_next(ex->send_map, &pos, &proc, &c))
  {
    // Move the sent elements to the back.
    int stride = msg->stride;
    if (msg->type == MPI_REAL_T)
    {
      real_t* array = (real_t*)orig_buffer;
      for (int i = 0; i < c->num_indices; ++i)
      {
        for (int s = 0; s < stride; ++s)
        {
          array[stride*c->indices[i]+s] = array[*count-1-num_sent];
          ++num_sent;
        }
      }
    }
    else if (msg->type == MPI_DOUBLE)
    {
      double* array = (double*)orig_buffer;
      for (int i = 0; i < c->num_indices; ++i)
      {
        for (int s = 0; s < stride; ++s)
        {
          array[stride*c->indices[i]+s] = array[*count-1-num_sent];
          ++num_sent;
        }
      }
    }
    else 
    {
      // FIXME: What about other data types???
      polymec_error("exchanger_finish_transfer: unimplemented data type.");
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
  return serializer_new("exchanger", ex_size, ex_read, ex_write, DTOR(exchanger_free));
}

// This helper sets up the exchanger ex so that it can distribute data from the 
// root process (0) according to the global partition vector. The partition 
// vector is NULL on all nonzero ranks and defined on rank 0.
exchanger_t* create_distributor(MPI_Comm comm, 
                                int64_t* global_partition,
                                int num_global_vertices)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank == 0) || (global_partition == NULL));
  ASSERT((rank == 0) || (num_global_vertices == 0));

  exchanger_t* distributor = exchanger_new(comm);
  
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
      int vertices[num_vertices_to_send[p]], k = 0, p_tag = p;
      for (int i = 0; i < num_global_vertices; ++i)
      {
        if (global_partition[i] == p)
          vertices[k++] = i;
      }
      MPI_Send(vertices, num_vertices_to_send[p], MPI_INT, p, p_tag, comm);
      exchanger_set_send(distributor, p, vertices, num_vertices_to_send[p], true);
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
    int local_vertices[num_local_vertices], p_tag = rank;
    MPI_Status status;
    MPI_Recv(local_vertices, num_local_vertices, MPI_INT, 0, p_tag, comm, &status);

    // Now register all the vertices we're receiving with the migrator.
    ASSERT(num_local_vertices > 0);
    exchanger_set_receive(distributor, 0, local_vertices, num_local_vertices, true);
  }

  return distributor;
}

// This helper sets up the exchanger ex so that it can migrate data from the 
// current process according to the local partition vector.
exchanger_t* create_migrator(MPI_Comm comm,
                             int64_t* local_partition,
                             int num_vertices)
{
  exchanger_t* migrator = exchanger_new(comm);
  
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

  // Send and receive the actual vertices.
  int** receive_vertices = polymec_malloc(sizeof(int*) * nprocs);
  memset(receive_vertices, 0, sizeof(int*) * nprocs);
  int num_requests = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    if ((rank != p) && (num_vertices_to_receive[p] > 0)) ++num_requests;
    if ((rank != p) && (num_vertices_to_send[p] > 0)) ++num_requests;
  }
  if (num_requests > 0)
  {
    MPI_Request requests[num_requests];
    int r = 0;
    for (int p = 0; p < nprocs; ++p)
    {
      if (num_vertices_to_receive[p] > 0)
      {
        receive_vertices[p] = polymec_malloc(sizeof(int) * num_vertices_to_receive[p]);
        if (p != rank)
        {
          // Note that we use this rank as our tag.
          MPI_Irecv(receive_vertices[p], num_vertices_to_receive[p], MPI_INT, p, rank, comm, &requests[r++]);
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
        {
          // Note that we use the destination rank p as our tag.
          exchanger_set_send(migrator, p, send_vertices, num_vertices_to_send[p], true);
          MPI_Isend(send_vertices, num_vertices_to_send[p], MPI_INT, p, p, comm, &requests[r++]);
        }
        else
          memcpy(receive_vertices[rank], send_vertices, sizeof(int) * num_vertices_to_receive[rank]);
      }
    }
    ASSERT(r == num_requests);

    // Wait for exchanges to finish. We can't use MPI_Waitall here because 
    // not all processes necessary participate.
    MPI_Status statuses[num_requests];
    int finished[num_requests];
    memset(finished, 0, num_requests*sizeof(int));
    while (true)
    {
      bool all_finished = true;
      for (int i = 0; i < num_requests; ++i)
      {
        if (!finished[i])
        {
          if (MPI_Test(&requests[i], &finished[i], &statuses[i]) != MPI_SUCCESS)
            polymec_error("create_migrator: Internal error.");
          if (!finished[i]) all_finished = false;
        }
      }

      // If the transmissions have finished at this point, we 
      // can break out of the loop. 
      if (all_finished) break;
    } 

    // Now register all the vertices we're receiving with the migrator.
    for (int p = 0; p < nprocs; ++p)
    {
      if (num_vertices_to_receive[p] > 0)
      {
        ASSERT(receive_vertices[p] != NULL);
        if (rank != p)
          exchanger_set_receive(migrator, p, receive_vertices[p], num_vertices_to_receive[p], true);
        polymec_free(receive_vertices[p]);
      }
    }
    polymec_free(receive_vertices);
  }

  return migrator;
}

