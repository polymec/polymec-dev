// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/exchanger.h"
#include "core/timer.h"
#include "core/unordered_map.h"
#include "core/unordered_set.h"

// Here's a virtual table for exchanger reduction operators.
typedef struct
{
  double    (*reduce_double)(void* context, double* values, int* processes, size_t num_values);
  float     (*reduce_float)(void* context, float* values, int* processes, size_t num_values);
  double complex (*reduce_double_complex)(void* context, double complex* values, int* processes, size_t num_values);
  float complex (*reduce_float_complex)(void* context, float complex* values, int* processes, size_t num_values);
  int       (*reduce_int)(void* context, int* values, int* processes, size_t num_values);
  long      (*reduce_long)(void* context, long* values, int* processes, size_t num_values);
  long long (*reduce_long_long)(void* context, long long* values, int* processes, size_t num_values);
  uint64_t  (*reduce_uint64)(void* context, uint64_t* values, int* processes, size_t num_values);
  int64_t   (*reduce_int64)(void* context, int64_t* values, int* processes, size_t num_values);
  char      (*reduce_char)(void* context, char* values, int* processes, size_t num_values);
  uint8_t   (*reduce_byte)(void* context, uint8_t* values, int* processes, size_t num_values);
  void (*dtor)(void* context);
} exchanger_reducer_vtable;

// And here's the reduction operator itself.
struct exchanger_reducer_t
{
  char* name;
  void* context;
  exchanger_reducer_vtable vtable;
};

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

static size_t mpi_size(MPI_Datatype type)
{
  size_t size = 0;
  if (type == MPI_DOUBLE)
    size = sizeof(double);
  else if (type == MPI_FLOAT)
    size = sizeof(float);
  else if (type == MPI_COMPLEX_T)
    size = sizeof(complex_t);
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

#define PACK(msg, mpi_type, c_type) \
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

    PACK(msg, MPI_DOUBLE, double)
    else PACK(msg, MPI_FLOAT, float)
    else PACK(msg, MPI_INT, int)
    else PACK(msg, MPI_LONG, long)
    else PACK(msg, MPI_LONG_LONG, long long)
    else PACK(msg, MPI_UINT64_T, uint64_t)
    else PACK(msg, MPI_INT64_T, int64_t)
    else PACK(msg, MPI_CHAR, char)
    else PACK(msg, MPI_BYTE, uint8_t)
    else PACK(msg, MPI_C_DOUBLE_COMPLEX, double complex)
    else PACK(msg, MPI_C_FLOAT_COMPLEX, float complex)
    else polymec_error("mpi_message_pack: unsupported type!");
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
#undef PACK

#define UNPACK(msg, mpi_type, c_type) \
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
    UNPACK(msg, MPI_DOUBLE, double)
    else UNPACK(msg, MPI_FLOAT, float)
    else UNPACK(msg, MPI_INT, int)
    else UNPACK(msg, MPI_LONG, long)
    else UNPACK(msg, MPI_LONG_LONG, long long)
    else UNPACK(msg, MPI_UINT64_T, uint64_t)
    else UNPACK(msg, MPI_INT64_T, int64_t)
    else UNPACK(msg, MPI_CHAR, char)
    else UNPACK(msg, MPI_BYTE, uint8_t)
    else UNPACK(msg, MPI_C_DOUBLE_COMPLEX, double complex)
    else UNPACK(msg, MPI_C_FLOAT_COMPLEX, float complex)
    else polymec_error("mpi_message_unpack: unsupported type!");
    ++i;
  }
}

#define UNPACK_AND_REDUCE(msg, mpi_type, c_type, method) \
  if (msg->type == mpi_type) \
  { \
    c_type* src = msg->receive_buffers[i]; \
    c_type* dest = data; \
    for (int j = 0; j < msg->receive_buffer_sizes[i]; ++j) \
    { \
      int index = recv_indices[j]; \
      int_array_t** procs_p = agg_map_get(aggregates, index); \
      if (procs_p != NULL) \
      { \
        int_array_t* procs = *procs_p; \
        size_t num_values = procs->size; \
        c_type values[stride][num_values]; \
        for (size_t n = 0; n < num_values; ++n) \
        { \
          int proc1 = procs->data[n]; \
          if (proc == proc1) \
          { \
            for (int s = 0; s < stride; ++s) \
              values[s][n] = src[stride*j+s]; \
          } \
          else \
          { \
            int k = 0; \
            while (msg->source_procs[k] != proc1) ++k; \
            ASSERT(k < msg->num_receives); \
            c_type* src1 = msg->receive_buffers[k]; \
            exchanger_channel_t* ch = *exchanger_map_get(receive_map, proc1); \
            int l = 0; \
            while (ch->indices[l] != index) ++l; \
            ASSERT(l < ch->num_indices); \
            for (int s = 0; s < stride; ++s) \
              values[s][n] = src1[stride*l+s]; \
          } \
        } \
        for (int s = 0; s < stride; ++s) \
        { \
          dest[receive_offset+stride*index+s] = \
            reducer->vtable.method(reducer->context, values[s], procs->data, num_values); \
        }\
      } \
      else \
      { \
        for (int s = 0; s < stride; ++s) \
          dest[receive_offset+stride*index+s] = src[stride*j+s]; \
      } \
    } \
  }

DEFINE_UNORDERED_MAP(agg_map, int, int_array_t*, int_hash, int_equals)

static void mpi_message_unpack_and_reduce(mpi_message_t* msg, 
                                          void* data, 
                                          ssize_t receive_offset,
                                          exchanger_map_t* receive_map,
                                          exchanger_reducer_t* reducer,
                                          agg_map_t* aggregates)
{
  int pos = 0, proc, i = 0;
  exchanger_channel_t* c;
  while (exchanger_map_next(receive_map, &pos, &proc, &c))
  {
    int* recv_indices = c->indices;
    int stride = msg->stride;
    UNPACK_AND_REDUCE(msg, MPI_DOUBLE, double, reduce_double)
    else UNPACK_AND_REDUCE(msg, MPI_FLOAT, float, reduce_float)
    else UNPACK_AND_REDUCE(msg, MPI_INT, int, reduce_int)
    else UNPACK_AND_REDUCE(msg, MPI_LONG, long, reduce_long)
    else UNPACK_AND_REDUCE(msg, MPI_LONG_LONG, long long, reduce_long_long)
    else UNPACK_AND_REDUCE(msg, MPI_UINT64_T, uint64_t, reduce_uint64)
    else UNPACK_AND_REDUCE(msg, MPI_INT64_T, int64_t, reduce_int64)
    else UNPACK_AND_REDUCE(msg, MPI_CHAR, char, reduce_char)
    else UNPACK_AND_REDUCE(msg, MPI_BYTE, uint8_t, reduce_byte)
    else UNPACK_AND_REDUCE(msg, MPI_C_DOUBLE_COMPLEX, double complex, reduce_double_complex)
    else UNPACK_AND_REDUCE(msg, MPI_C_FLOAT_COMPLEX, float complex, reduce_float_complex)
    else polymec_error("mpi_message_unpack_and_reduce: unsupported type!");
    ++i;
  }
}

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
  else if (msg->type == MPI_COMPLEX_T)
    strcpy(typeStr, "complex");
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

  // Reducer function and aggregates map linking receive index -> procs.
  exchanger_reducer_t* reducer;
  agg_map_t* agg_procs;
};

static void exchanger_clear(exchanger_t* ex)
{
  agg_map_clear(ex->agg_procs);
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
  agg_map_free(ex->agg_procs);
}

static void init_reducers(void);
exchanger_t* exchanger_new_with_rank(MPI_Comm comm, int rank)
{
  // If we haven't yet done so, initialize our built-in reducers.
  static bool reducers_initialized = false;
  if (!reducers_initialized)
  {
    init_reducers();
    reducers_initialized = true;
  }

  exchanger_t* ex = polymec_refcounted_malloc(sizeof(exchanger_t), exchanger_free);
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
  ex->pending_msgs = polymec_calloc(ex->pending_msg_cap, sizeof(mpi_message_t*));
  ex->orig_buffers = polymec_calloc(ex->pending_msg_cap, sizeof(void*));
  ex->transfer_counts = polymec_calloc(ex->pending_msg_cap, sizeof(int*));
  ex->max_send = -1;
  ex->max_receive = -1;
  ex->reducer = NULL;
  ex->agg_procs = agg_map_new();

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

void exchanger_set_send(exchanger_t* ex, 
                        int remote_process, 
                        int* indices, 
                        int num_indices, 
                        bool copy_indices)
{
  ASSERT(remote_process >= 0);
  ASSERT(remote_process < ex->nprocs);
  ASSERT(num_indices >= 0);

  if (num_indices > 0)
  {
    exchanger_channel_t* c = exchanger_channel_new(num_indices, indices, copy_indices);
    exchanger_map_insert_with_kv_dtor(ex->send_map, remote_process, c, delete_map_entry);

    if (remote_process > ex->max_send)
      ex->max_send = remote_process;
  }
}

void exchanger_set_sends(exchanger_t* ex, exchanger_proc_map_t* send_map)
{
  int pos = 0, send_proc;
  int_array_t* send_indices;
  while (exchanger_proc_map_next(send_map, &pos, &send_proc, &send_indices))
    exchanger_set_send(ex, send_proc, send_indices->data, (int)send_indices->size, true);
  exchanger_proc_map_free(send_map);
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

static void set_receive(exchanger_t* ex, 
                        int remote_process, 
                        int* indices, 
                        int num_indices, 
                        bool copy_indices)
{
  if (num_indices > 0)
  {
    // Set up the mapping in our channel.
    exchanger_channel_t* c = exchanger_channel_new(num_indices, indices, copy_indices);
    exchanger_map_insert_with_kv_dtor(ex->receive_map, remote_process, c, delete_map_entry);

    // Update the maximum remote process if needed.
    if (remote_process > ex->max_receive)
      ex->max_receive = remote_process;
  }
}

static void find_aggregates(exchanger_t* ex)
{
  START_FUNCTION_TIMER();

  // Clear existing aggregate entries.
  agg_map_clear(ex->agg_procs);

  // Keep track of receive indices we've already encountered here.
  int_int_unordered_map_t* receive_indices = int_int_unordered_map_new();

  // Scour the receive map for aggregated values.
  int pos = 0, proc;
  exchanger_channel_t* c;
  while (exchanger_map_next(ex->receive_map, &pos, &proc, &c))
  {
    for (int i = 0; i < c->num_indices; ++i)
    {
      int index = c->indices[i];
      int* proc_p = int_int_unordered_map_get(receive_indices, index);
      if (proc_p != NULL)
      {
        // We've already encountered this receive index, so start 
        // aggregating the processes.
        int_array_t** procs_p = agg_map_get(ex->agg_procs, index);
        int_array_t* procs = NULL;
        if (procs_p == NULL)
        {
          procs = int_array_new();
          agg_map_insert_with_v_dtor(ex->agg_procs, index, procs, int_array_free);
          int_array_append(procs, *proc_p);
        }
        else
          procs = *procs_p;
        int_array_append(procs, proc);
      }
      else
        int_int_unordered_map_insert(receive_indices, index, proc);
    }
  }
  int_int_unordered_map_free(receive_indices);
  STOP_FUNCTION_TIMER();
}

void exchanger_set_receive(exchanger_t* ex, 
                           int remote_process, 
                           int* indices, 
                           int num_indices, 
                           bool copy_indices)
{
  ASSERT(remote_process >= 0);
  ASSERT(remote_process < ex->nprocs);
  ASSERT(num_indices >= 0);

  if (num_indices > 0)
  {
    set_receive(ex, remote_process, indices, num_indices, copy_indices);
    find_aggregates(ex);
  }
}

void exchanger_set_receives(exchanger_t* ex, exchanger_proc_map_t* recv_map)
{
  int pos = 0, recv_proc;
  int_array_t* recv_indices;
  while (exchanger_proc_map_next(recv_map, &pos, &recv_proc, &recv_indices))
    set_receive(ex, recv_proc, recv_indices->data, (int)recv_indices->size, true);   
  find_aggregates(ex);
  exchanger_proc_map_free(recv_map);
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
  log_debug("exchanger_verify: Checking connectivity.");

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
        STOP_FUNCTION_TIMER();
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
        STOP_FUNCTION_TIMER();
        return false;
      }
    }
  }
  polymec_free(neighbors_for_proc);
  log_debug("exchanger_verify: Connectivity verified successfully.");
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

static int exchanger_send_message(exchanger_t* ex, mpi_message_t* msg)
{
  START_FUNCTION_TIMER();

  int j = 0, idest_local = -1;
  for (int i = 0; i < msg->num_receives; ++i)
  {
#if POLYMEC_HAVE_MPI
    if (ex->rank != msg->source_procs[i])
    {
      // If we are expecting data, post an asynchronous receive. 
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
    else // Mark down the index of the local copy.
#endif
      idest_local = i;
  }

  // Send the data asynchronously. 
  for (int i = 0; i < msg->num_sends; ++i)
  {
#if POLYMEC_HAVE_MPI
    if (ex->rank != msg->dest_procs[i])
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
    else 
#endif
    {
      // Do any local copying.
      ASSERT(idest_local != -1); // no local destination for this source?
      ASSERT(msg->receive_buffer_sizes[idest_local] == msg->send_buffer_sizes[i]);
      memcpy(msg->receive_buffers[idest_local], msg->send_buffers[i], 
             mpi_size(msg->type) * msg->stride * msg->send_buffer_sizes[i]);
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
  STOP_FUNCTION_TIMER();
  return token;
}

int exchanger_start_exchange(exchanger_t* ex, void* data, int stride, int tag, MPI_Datatype type)
{
  START_FUNCTION_TIMER();

  // If we aggregate any data, we'd better have a reducer.
  if ((ex->agg_procs->size > 0) && (ex->reducer == NULL))
  {
    polymec_error("exchanger_start_exchange: exchanger aggregates data, "
                  "but has no reducer set.");
  }

  // Create a message for this array.
  mpi_message_t* msg = mpi_message_new(type, stride, tag);
  mpi_message_pack(msg, data, ex->send_offset, ex->send_map, ex->receive_map);
  
  // Begin the transmission and allocate a token for it.
  int token = exchanger_send_message(ex, msg);
  ex->orig_buffers[token] = data;
  STOP_FUNCTION_TIMER();
  return token;
}

static int exchanger_waitall(exchanger_t* ex, mpi_message_t* msg)
{
#if POLYMEC_HAVE_MPI
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
#else
  return 0;
#endif
}

void exchanger_finish_exchange(exchanger_t* ex, int token) 
{
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT(token < ex->num_pending_msgs);
  ASSERT(ex->transfer_counts[token] == NULL); // We can't finish a transfer as an exchange!

  // Retrieve the message for the given token.
  mpi_message_t* msg = ex->pending_msgs[token];
  int err = exchanger_waitall(ex, msg);
  if (err == MPI_SUCCESS) 
  {
    void* orig_buffer = ex->orig_buffers[token];
    if (ex->agg_procs->size == 0)
    {
      // Copy the received data into the original array.
      mpi_message_unpack(msg, orig_buffer, ex->receive_offset, ex->receive_map);
    }
    else
    {
      ASSERT(ex->reducer != NULL);
      // Copy and reduce the received data into the original array.
      mpi_message_unpack_and_reduce(msg, orig_buffer, ex->receive_offset, 
                                    ex->receive_map, ex->reducer,
                                    ex->agg_procs);
    }
  }

  // Pull the message out of our list of pending messages and delete it.
  ex->pending_msgs[token] = NULL;
  ex->orig_buffers[token] = NULL;
  mpi_message_free(msg);
  STOP_FUNCTION_TIMER();
}

void exchanger_fprintf(exchanger_t* ex, FILE* stream)
{
  fprintf(stream, "Exchanger(");
  if (ex->comm == MPI_COMM_WORLD)
    fprintf(stream, "MPI_COMM_WORLD");
  else if (ex->comm == MPI_COMM_SELF)
    fprintf(stream, "MPI_COMM_SELF");
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
}

static size_t ex_size(void* obj)
{
  exchanger_t* ex = obj;
  size_t size = 2 * sizeof(int);
  int pos = 0, proc, *indices, num_indices;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
    size += (2 + num_indices) * sizeof(int);
  pos = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
    size += (2 + num_indices) * sizeof(int);
  return size;
}

static void* ex_read(byte_array_t* bytes, size_t* offset)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  exchanger_t* ex = exchanger_new(comm);
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
  return ex;
}

static void ex_write(void* obj, byte_array_t* bytes, size_t* offset)
{
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
}

serializer_t* exchanger_serializer()
{
  return serializer_new("exchanger", ex_size, ex_read, ex_write, NULL);
}

// This helper function gets rid of a lot of boilerplate code.
void exchanger_proc_map_add_index(exchanger_proc_map_t* map, int process, int index)
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
  int_array_append(indices, index);
}

static size_t epm_size(void* obj)
{
  exchanger_proc_map_t* epm = obj;
  size_t size = sizeof(int);
  int pos = 0, proc;
  int_array_t* indices;
  while (exchanger_proc_map_next(epm, &pos, &proc, &indices))
    size += (2 + indices->size) * sizeof(int);
  return size;
}

static void* epm_read(byte_array_t* bytes, size_t* offset)
{
  exchanger_proc_map_t* epm = exchanger_proc_map_new();
  int size;
  byte_array_read_ints(bytes, 1, &size, offset);
  for (int i = 0; i < size; ++i)
  {
    int proc, num_indices;
    byte_array_read_ints(bytes, 1, &proc, offset);
    byte_array_read_ints(bytes, 1, &num_indices, offset);
    int indices[num_indices];
    byte_array_read_ints(bytes, num_indices, indices, offset);
    for (int j = 0; j < num_indices; ++j)
      exchanger_proc_map_add_index(epm, proc, indices[j]);
  }
  return epm;
}

static void epm_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  exchanger_proc_map_t* epm = obj;

  int pos = 0, proc;
  int_array_t* indices;
  int size = (int)epm->size;
  byte_array_write_ints(bytes, 1, &size, offset);
  while (exchanger_proc_map_next(epm, &pos, &proc, &indices))
  {
    byte_array_write_ints(bytes, 1, &proc, offset);
    int num_indices = (int)indices->size;
    byte_array_write_ints(bytes, 1, &num_indices, offset);
    byte_array_write_ints(bytes, num_indices, indices->data, offset);
  }
}

serializer_t* exchanger_proc_map_serializer()
{
  return serializer_new("exchanger_proc_map", epm_size, epm_read, epm_write, NULL);
}

bool exchanger_aggregates_data(exchanger_t* ex)
{
  return (ex->agg_procs->size > 0);
}

void exchanger_set_reducer(exchanger_t* ex,
                           exchanger_reducer_t* reducer)
{
  ex->reducer = reducer;
}

//------------------------------------------------------------------------
//                         Built-in reducers
//------------------------------------------------------------------------

static exchanger_reducer_t* exchanger_reducer_new(const char* name,
                                                  void* context,
                                                  exchanger_reducer_vtable vtable)
{
  exchanger_reducer_t* reducer = polymec_malloc(sizeof(exchanger_reducer_t));
  reducer->name = string_dup(name);
  reducer->context = context;
  reducer->vtable = vtable;
  return reducer;
}

static void exchanger_reducer_free(exchanger_reducer_t* reducer)
{
  if ((reducer->context != NULL) && (reducer->vtable.dtor != NULL))
    reducer->vtable.dtor(reducer->context);
  string_free(reducer->name);
  polymec_free(reducer);
}

// Built-in reducers.
exchanger_reducer_t* EXCHANGER_SUM = NULL;
exchanger_reducer_t* EXCHANGER_PRODUCT = NULL;
exchanger_reducer_t* EXCHANGER_MIN = NULL;
exchanger_reducer_t* EXCHANGER_MAX = NULL;
exchanger_reducer_t* EXCHANGER_MIN_RANK = NULL;
exchanger_reducer_t* EXCHANGER_MAX_RANK = NULL;

static void free_reducers(void)
{
  exchanger_reducer_free(EXCHANGER_SUM);
  exchanger_reducer_free(EXCHANGER_PRODUCT);
  exchanger_reducer_free(EXCHANGER_MIN);
  exchanger_reducer_free(EXCHANGER_MAX);
  exchanger_reducer_free(EXCHANGER_MIN_RANK);
  exchanger_reducer_free(EXCHANGER_MAX_RANK);
}

#define DEFINE_SUM(func_name, c_type, zero) \
static c_type func_name(void* context, c_type* values, int* processes, size_t num_values) \
{ \
  c_type result = zero; \
  for (size_t i = 0; i < num_values; ++i) \
    result += values[i]; \
  return result; \
} 
DEFINE_SUM(sum_double, double, 0.0)
DEFINE_SUM(sum_float, float, 0.0)
DEFINE_SUM(sum_double_complex, double complex, CMPLX(0.0, 0.0))
DEFINE_SUM(sum_float_complex, float complex, CMPLXF(0.0, 0.0))
DEFINE_SUM(sum_int, int, 0)
DEFINE_SUM(sum_long, long, 0)
DEFINE_SUM(sum_long_long, long long, 0)
DEFINE_SUM(sum_uint64, uint64_t, 0)
DEFINE_SUM(sum_int64, int64_t, 0)
DEFINE_SUM(sum_char, char, 0)
DEFINE_SUM(sum_byte, uint8_t, 0)

#define DEFINE_PRODUCT(func_name, c_type, one) \
static c_type func_name(void* context, c_type* values, int* processes, size_t num_values) \
{ \
  c_type result = one; \
  for (size_t i = 0; i < num_values; ++i) \
    result *= values[i]; \
  return result; \
} 
DEFINE_PRODUCT(prod_double, double, 1.0)
DEFINE_PRODUCT(prod_float, float, 1.0)
DEFINE_PRODUCT(prod_double_complex, double complex, CMPLX(1.0, 0.0))
DEFINE_PRODUCT(prod_float_complex, float complex, CMPLXF(1.0, 0.0))
DEFINE_PRODUCT(prod_int, int, 1)
DEFINE_PRODUCT(prod_long, long, 1)
DEFINE_PRODUCT(prod_long_long, long long, 1)
DEFINE_PRODUCT(prod_uint64, uint64_t, 1)
DEFINE_PRODUCT(prod_int64, int64_t, 1)
DEFINE_PRODUCT(prod_char, char, 1)
DEFINE_PRODUCT(prod_byte, uint8_t, 1)

#define DEFINE_MINMAX(func_name, c_type, start, min_or_max) \
static c_type func_name(void* context, c_type* values, int* processes, size_t num_values) \
{ \
  c_type result = start; \
  for (size_t i = 0; i < num_values; ++i) \
    result = min_or_max(result, values[i]); \
  return result; \
} 
DEFINE_MINMAX(min_double, double, DBL_MAX, MIN)
DEFINE_MINMAX(min_float, float, FLT_MAX, MIN)
DEFINE_MINMAX(min_int, int, INT_MAX, MIN)
DEFINE_MINMAX(min_long, long, INT_MAX, MIN)
DEFINE_MINMAX(min_long_long, long long, INT_MAX, MIN)
DEFINE_MINMAX(min_uint64, uint64_t, INT_MAX, MIN)
DEFINE_MINMAX(min_int64, int64_t, INT_MAX, MIN)
DEFINE_MINMAX(min_char, char, 127, MIN)
DEFINE_MINMAX(min_byte, uint8_t, 255, MIN)

DEFINE_MINMAX(max_double, double, -DBL_MAX, MAX)
DEFINE_MINMAX(max_float, float, -FLT_MAX, MAX)
DEFINE_MINMAX(max_int, int, -INT_MAX, MAX)
DEFINE_MINMAX(max_long, long, -INT_MAX, MAX)
DEFINE_MINMAX(max_long_long, long long, -INT_MAX, MAX)
DEFINE_MINMAX(max_uint64, uint64_t, 0, MAX)
DEFINE_MINMAX(max_int64, int64_t, -INT_MAX, MAX)
DEFINE_MINMAX(max_char, char, -127, MAX)
DEFINE_MINMAX(max_byte, uint8_t, 0, MAX)

// For complex numbers, the MIN and MAX functions use the modulus.
static double complex min_double_complex(void* context, double complex* values, int* processes, size_t num_values)
{
  double complex result = CMPLX(DBL_MAX, 0.0); 
  for (size_t i = 0; i < num_values; ++i) 
    result = MIN(cabs(result), cabs(values[i])); 
  return result; 
}

static double complex max_double_complex(void* context, double complex* values, int* processes, size_t num_values)
{
  double complex result = CMPLX(0.0, 0.0); 
  for (size_t i = 0; i < num_values; ++i) 
    result = MAX(cabs(result), cabs(values[i])); 
  return result; 
}

static float complex min_float_complex(void* context, float complex* values, int* processes, size_t num_values)
{
  float complex result = CMPLXF(FLT_MAX, 0.0); 
  for (size_t i = 0; i < num_values; ++i) 
    result = MIN(cabsf(result), cabsf(values[i])); 
  return result; 
}

static float complex max_float_complex(void* context, float complex* values, int* processes, size_t num_values)
{
  float complex result = CMPLXF(0.0, 0.0); 
  for (size_t i = 0; i < num_values; ++i) 
    result = MAX(cabsf(result), cabsf(values[i])); 
  return result; 
}

#define DEFINE_MIN_RANK(func_name, c_type) \
static c_type func_name(void* context, c_type* values, int* processes, size_t num_values) \
{ \
  int proc = processes[0]; \
  c_type result = values[0]; \
  for (size_t i = 1; i < num_values; ++i) \
  { \
    if (processes[i] < proc) \
    { \
      proc = processes[i]; \
      result = values[i]; \
    } \
  } \
  return result; \
} 
DEFINE_MIN_RANK(minr_double, double)
DEFINE_MIN_RANK(minr_float, float)
DEFINE_MIN_RANK(minr_double_complex, double complex)
DEFINE_MIN_RANK(minr_float_complex, float complex)
DEFINE_MIN_RANK(minr_int, int)
DEFINE_MIN_RANK(minr_long, long)
DEFINE_MIN_RANK(minr_long_long, long long)
DEFINE_MIN_RANK(minr_uint64, uint64_t)
DEFINE_MIN_RANK(minr_int64, int64_t)
DEFINE_MIN_RANK(minr_char, char)
DEFINE_MIN_RANK(minr_byte, uint8_t)

#define DEFINE_MAX_RANK(func_name, c_type) \
static c_type func_name(void* context, c_type* values, int* processes, size_t num_values) \
{ \
  int proc = processes[0]; \
  c_type result = values[0]; \
  for (size_t i = 1; i < num_values; ++i) \
  { \
    if (processes[i] > proc) \
    { \
      proc = processes[i]; \
      result = values[i]; \
    } \
  } \
  return result; \
} 
DEFINE_MAX_RANK(maxr_double, double)
DEFINE_MAX_RANK(maxr_float, float)
DEFINE_MAX_RANK(maxr_double_complex, double complex)
DEFINE_MAX_RANK(maxr_float_complex, float complex)
DEFINE_MAX_RANK(maxr_int, int)
DEFINE_MAX_RANK(maxr_long, long)
DEFINE_MAX_RANK(maxr_long_long, long long)
DEFINE_MAX_RANK(maxr_uint64, uint64_t)
DEFINE_MAX_RANK(maxr_int64, int64_t)
DEFINE_MAX_RANK(maxr_char, char)
DEFINE_MAX_RANK(maxr_byte, uint8_t)

static void init_reducers(void)
{
  exchanger_reducer_vtable sum_vtable = {.reduce_double = sum_double,
                                         .reduce_float = sum_float,
                                         .reduce_double_complex = sum_double_complex,
                                         .reduce_float_complex = sum_float_complex,
                                         .reduce_int = sum_int,
                                         .reduce_long = sum_long,
                                         .reduce_long_long = sum_long_long,
                                         .reduce_uint64 = sum_uint64,
                                         .reduce_int64 = sum_int64,
                                         .reduce_char = sum_char,
                                         .reduce_byte = sum_byte};
  EXCHANGER_SUM = exchanger_reducer_new("EXCHANGER_SUM", NULL, sum_vtable);

  exchanger_reducer_vtable prod_vtable = {.reduce_double = prod_double,
                                          .reduce_float = prod_float,
                                          .reduce_double_complex = prod_double_complex,
                                          .reduce_float_complex = prod_float_complex,
                                          .reduce_int = prod_int,
                                          .reduce_long = prod_long,
                                          .reduce_long_long = prod_long_long,
                                          .reduce_uint64 = prod_uint64,
                                          .reduce_int64 = prod_int64,
                                          .reduce_char = prod_char,
                                          .reduce_byte = prod_byte};
  EXCHANGER_PRODUCT = exchanger_reducer_new("EXCHANGER_PRODUCT", NULL, prod_vtable);

  exchanger_reducer_vtable min_vtable = {.reduce_double = min_double,
                                         .reduce_float = min_float,
                                         .reduce_double_complex = min_double_complex,
                                         .reduce_float_complex = min_float_complex,
                                         .reduce_int = min_int,
                                         .reduce_long = min_long,
                                         .reduce_long_long = min_long_long,
                                         .reduce_uint64 = min_uint64,
                                         .reduce_int64 = min_int64,
                                         .reduce_char = min_char,
                                         .reduce_byte = min_byte};
  EXCHANGER_MIN = exchanger_reducer_new("EXCHANGER_MIN", NULL, min_vtable);

  exchanger_reducer_vtable max_vtable = {.reduce_double = max_double,
                                         .reduce_float = max_float,
                                         .reduce_double_complex = max_double_complex,
                                         .reduce_float_complex = max_float_complex,
                                         .reduce_int = max_int,
                                         .reduce_long = max_long,
                                         .reduce_long_long = max_long_long,
                                         .reduce_uint64 = max_uint64,
                                         .reduce_int64 = max_int64,
                                         .reduce_char = max_char,
                                         .reduce_byte = max_byte};
  EXCHANGER_MAX = exchanger_reducer_new("EXCHANGER_MAX", NULL, max_vtable);

  exchanger_reducer_vtable minr_vtable = {.reduce_double = minr_double,
                                          .reduce_float = minr_float,
                                          .reduce_double_complex = minr_double_complex,
                                          .reduce_float_complex = minr_float_complex,
                                          .reduce_int = minr_int,
                                          .reduce_long = minr_long,
                                          .reduce_long_long = minr_long_long,
                                          .reduce_uint64 = minr_uint64,
                                          .reduce_int64 = minr_int64,
                                          .reduce_char = minr_char,
                                          .reduce_byte = minr_byte};
  EXCHANGER_MIN_RANK = exchanger_reducer_new("EXCHANGER_MIN_RANK", NULL, minr_vtable);

  exchanger_reducer_vtable maxr_vtable = {.reduce_double = maxr_double,
                                          .reduce_float = maxr_float,
                                          .reduce_double_complex = maxr_double_complex,
                                          .reduce_float_complex = maxr_float_complex,
                                          .reduce_int = maxr_int,
                                          .reduce_long = maxr_long,
                                          .reduce_long_long = maxr_long_long,
                                          .reduce_uint64 = maxr_uint64,
                                          .reduce_int64 = maxr_int64,
                                          .reduce_char = maxr_char,
                                          .reduce_byte = maxr_byte};
  EXCHANGER_MAX_RANK = exchanger_reducer_new("EXCHANGER_MAX_RANK", NULL, maxr_vtable);

  polymec_atexit(free_reducers);
}

