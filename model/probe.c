// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "model/probe.h"
#include "model/model.h"

probe_data_t* probe_data_new(int rank, size_t* shape)
{
  ASSERT(rank >= 0);
  probe_data_t* data = polymec_malloc(sizeof(probe_data_t));
  data->rank = rank;
  data->shape = polymec_malloc(sizeof(size_t) * rank);
  size_t n = 1;
  for (int i = 0; i < rank; ++i)
  {
    data->shape[i] = shape[i];
    n *= shape[i];
  }
  data->data = polymec_malloc(sizeof(real_t) * n);
  return data;
}

void probe_data_free(probe_data_t* data)
{
  polymec_free(data->data);
  polymec_free(data->shape);
  polymec_free(data);
}

// Here's an acquisition callback.
typedef struct
{
  void* context;
  void (*function)(void* context, real_t t, probe_data_t* data);
  void (*dtor)(void* context);
} acq_callback_t;

static void free_callback_context(acq_callback_t callback)
{
  if ((callback.context != NULL) && (callback.dtor != NULL))
    callback.dtor(callback.context);
}

DEFINE_ARRAY(acq_callback_array, acq_callback_t)

struct probe_t
{
  char* name;
  char* data_name;
  int rank;
  size_t* shape;
  void* context;
  probe_vtable vtable;
  acq_callback_array_t* callbacks;
};

probe_t* probe_new(const char* name, 
                   const char* data_name,
                   int rank,
                   size_t* shape,
                   void* context, 
                   probe_vtable vtable)
{
  ASSERT(rank >= 0);
  ASSERT((shape != NULL) || (rank == 0));
  ASSERT(vtable.acquire != NULL);

  probe_t* probe = polymec_malloc(sizeof(probe_t));
  probe->name = string_dup(name);
  probe->data_name = string_dup(data_name);
  probe->rank = rank;
  probe->shape = polymec_malloc(sizeof(size_t) * rank);
  probe->context = context;
  for (int i = 0; i < rank; ++i)
  {
    ASSERT(shape[i] > 0);
    probe->shape[i] = shape[i];
  }
  probe->vtable = vtable;
  probe->callbacks = acq_callback_array_new();
  return probe;
}

void probe_free(probe_t* probe)
{
  acq_callback_array_free(probe->callbacks);
  string_free(probe->name);
  string_free(probe->data_name);
  if ((probe->context != NULL) && (probe->vtable.dtor != NULL))
    probe->vtable.dtor(probe->context);
  polymec_free(probe->shape);
  polymec_free(probe);
}

char* probe_name(probe_t* probe)
{
  return probe->name;
}

char* probe_data_name(probe_t* probe)
{
  return probe->data_name;
}

void* probe_context(probe_t* probe)
{
  return probe->context;
}

probe_data_t* probe_acquire(probe_t* probe, real_t t)
{
  // Allocate storage for and acquire the data.
  probe_data_t* data = probe_data_new(probe->rank, probe->shape);
  probe->vtable.acquire(probe->context, t, data);
  data->time = t;

  // Call any callbacks we have with the time and the data.
  for (size_t i = 0; i < probe->callbacks->size; ++i)
  {
    acq_callback_t* callback = &(probe->callbacks->data[i]);
    callback->function(callback->context, t, data);
  }

  // Return the data.
  return data;
}

void probe_set_model(probe_t* probe, model_t* model);
void probe_set_model(probe_t* probe, model_t* model)
{
  if (probe->vtable.set_model != NULL)
    probe->vtable.set_model(probe_context(probe), model_context(model));
}

void probe_postprocess(probe_t* probe, real_array_t* times, probe_data_array_t* data)
{
  if ((probe->vtable.postprocess != NULL) && (data != NULL))
    probe->vtable.postprocess(probe->context, times, data);
}

void probe_on_acquire(probe_t* probe, 
                      void* context, 
                      void (*function)(void* context, real_t t, probe_data_t* data),
                      void (*dtor)(void* context))
{
  acq_callback_t callback = {.context = context, .function = function, .dtor = dtor};
  acq_callback_array_append_with_dtor(probe->callbacks, callback, free_callback_context);
}

//------------------------------------------------------------------------
//                              Streaming nonsense
//------------------------------------------------------------------------

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <unistd.h>

typedef struct
{
  int socket_fd;
  struct sockaddr_in inet_addr;
  struct sockaddr_un unix_addr;
  int port;

  char* data_name;
  size_t data_size;
} stream_context_t;

static void free_stream(void* context)
{
  stream_context_t* stream = context;
  close(stream->socket_fd);
  polymec_free(stream);
}

static void stream_on_acquire(void* context, real_t t, probe_data_t* data)
{
  // The datagram consists of (in a sequence of bytes):
  // 1. The length of the probe's data name
  // 2. The characters in the probe's data name (including '\0')
  // 3. The time of the acquisition.
  // 4. The number of real numbers in the data.
  // 5. The data.
  stream_context_t* stream = context;
  size_t data_name_len = strlen(stream->data_name);
  size_t buff_size = sizeof(size_t) + sizeof(char) * (data_name_len+1) + sizeof(size_t) + sizeof(real_t) * (1 + stream->data_size);
  char buff[buff_size];
  memcpy(buff, &data_name_len, sizeof(size_t));
  size_t offset = sizeof(size_t);
  memcpy(&buff[offset], stream->data_name, sizeof(char) * data_name_len);
  offset += sizeof(char) * data_name_len;
  buff[offset] = '\0';
  offset += sizeof(char);
  memcpy(&buff[offset], &t, sizeof(real_t));
  offset += sizeof(real_t);
  memcpy(&buff[offset], &(stream->data_size), sizeof(size_t));
  offset += sizeof(size_t);
  memcpy(&buff[offset], data->data, sizeof(real_t) * stream->data_size);
  offset += sizeof(real_t) * stream->data_size;
  write(stream->socket_fd, buff, buff_size);
}

bool probe_stream_on_acquire(probe_t* probe, const char* destination, int port)
{
  stream_context_t* context = polymec_malloc(sizeof(stream_context_t));
  context->data_name = probe->data_name;
  context->data_size = 1;
  for (int i = 0; i < probe->rank; ++i)
    context->data_size *= probe->shape[i];

  // Try connecting via UDP first.
  int result;
  context->socket_fd = socket(AF_INET, SOCK_DGRAM, 0);
  if (context->socket_fd <= 0)
  {
    polymec_free(context);
    return false;
  }
  memset(&(context->inet_addr), 0, sizeof(struct sockaddr_in));
  context->inet_addr.sin_family = AF_INET;
  memcpy(&(context->inet_addr.sin_addr), destination, strlen(destination));
  context->inet_addr.sin_port = htons(port);
  result = connect(context->socket_fd, (struct sockaddr*)&context->inet_addr, sizeof(struct sockaddr_in));

  // If that didn't work, try a UNIX domain socket.
  if (result == -1)
  {
    context->socket_fd = socket(AF_INET, SOCK_DGRAM, 0);
    if (context->socket_fd <= 0)
    {
      polymec_free(context);
      return false;
    }
    memset(&(context->unix_addr), 0, sizeof(struct sockaddr_un));
    context->unix_addr.sun_family = AF_UNIX;
    memcpy(context->unix_addr.sun_path, destination, strlen(destination));
    result = connect(context->socket_fd, (struct sockaddr*)&context->unix_addr, sizeof(struct sockaddr_un));
  }

  if (result != -1)
  {
    probe_on_acquire(probe, context, stream_on_acquire, free_stream);
    return true;
  }
  else
  {
    polymec_free(context);
    return false;
  }
}

