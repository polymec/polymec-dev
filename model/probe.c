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
#include <errno.h>

typedef struct
{
  int socket_fd;
  struct sockaddr_in inet_addr;
  struct sockaddr_un unix_addr;

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
  static const char* header = "polymec-probe-stream";
  static size_t header_len = 20;

  stream_context_t* stream = context;
  size_t data_name_len = strlen(stream->data_name);
  size_t buff_size = sizeof(char)*header_len + sizeof(char) + 
                     sizeof(size_t) + sizeof(char) * data_name_len + sizeof(char) + 
                     sizeof(size_t) + sizeof(real_t) * (1 + stream->data_size) + sizeof(char);

  // Assemble the buffer.
  char buff[buff_size];

  // 1. The string "polymec-probe-stream" (excluding '\0').
  memcpy(buff, header, sizeof(char) * header_len); 
  size_t offset = header_len;

  // 2. A newline.
  buff[offset] = '\n';
  ++offset;

  // 3. A size_t containing the length of the probe's data name.
  memcpy(&buff[offset], &data_name_len, sizeof(size_t));
  offset += sizeof(size_t);

  // 4. The characters in the probe's data name (excluding '\0').
  memcpy(&buff[offset], stream->data_name, sizeof(char) * data_name_len);
  offset += sizeof(char) * data_name_len;

  // 5. A newline.
  buff[offset] = '\n';
  ++offset;

  // 6. A real_t containing the time of the acquisition.
  memcpy(&buff[offset], &t, sizeof(real_t));
  offset += sizeof(real_t);

  // 7. A size_t containing the number of real numbers in the data.
  memcpy(&buff[offset], &(stream->data_size), sizeof(size_t));
  offset += sizeof(size_t);

  // 8. A number of real_t data.
  memcpy(&buff[offset], data->data, sizeof(real_t) * stream->data_size);
  offset += sizeof(real_t) * stream->data_size;

  // 9. A newline.
  buff[offset] = '\n';
  ++offset;
  ASSERT(offset == buff_size);

  // Write the buffer to the stream.
  write(stream->socket_fd, buff, buff_size);
}

bool probe_stream_on_acquire(probe_t* probe, const char* destination, int port)
{
  stream_context_t* context = polymec_malloc(sizeof(stream_context_t));
  context->data_name = probe->data_name;
  context->data_size = 1;
  for (int i = 0; i < probe->rank; ++i)
    context->data_size *= probe->shape[i];

  // If we're given a valid port, try connecting via UDP first.
  char err_msg[129];
  if (port > 0)
  {
    context->socket_fd = socket(AF_INET, SOCK_DGRAM, 0);
    if (context->socket_fd > 0)
    {
      memset(&(context->inet_addr), 0, sizeof(context->inet_addr));
      context->inet_addr.sin_family = AF_INET;
      memcpy(&(context->inet_addr.sin_addr), destination, sizeof(char) * strlen(destination));
      context->inet_addr.sin_port = htons(port);
      int result = connect(context->socket_fd, (struct sockaddr*)&(context->inet_addr), sizeof(context->inet_addr));
      if (result != -1)
      {
        log_info("Probe %s: streaming %s to %s:%d via UDP", probe->name, probe->data_name, destination, port);
        probe_on_acquire(probe, context, stream_on_acquire, free_stream);
        return true;
      }
      else
      {
        // Clean up our mess.
        snprintf(err_msg, 128, "connect() returned errno %d", errno);
        close(context->socket_fd);
      }
    }
    else
      snprintf(err_msg, 128, "Failed to create an AF_INET socket.");
  }

  // If we're still here, UDP didn't work out, or we aren't doing UDP.
  // Try to use a UNIX socket.
  context->socket_fd = socket(AF_UNIX, SOCK_DGRAM, 0);
  if (context->socket_fd > 0)
  {
    memset(&(context->unix_addr), 0, sizeof(context->unix_addr));
    context->unix_addr.sun_family = AF_UNIX;
    memcpy(context->unix_addr.sun_path, destination, sizeof(char) * strlen(destination));
    int result = connect(context->socket_fd, (struct sockaddr*)&(context->unix_addr), sizeof(context->unix_addr));
    if (result != -1)
    {
      log_info("Probe %s: streaming %s to %s via UNIX sockets", probe->name, probe->data_name, destination);
      probe_on_acquire(probe, context, stream_on_acquire, free_stream);
      return true;
    }
    else
    {
      snprintf(err_msg, 128, "connect() returned errno %d", errno);
      close(context->socket_fd);
    }
  }
  else
    snprintf(err_msg, 128, "Failed to create an AF_UNIX socket.");

  // Better luck next time!
  log_info("Probe %s: Can't stream %s to %s (%s)", probe->name, probe->data_name, destination, err_msg);
  polymec_free(context);
  return false;
}

