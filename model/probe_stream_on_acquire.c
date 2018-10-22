// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/probe.h"
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>
#include <unistd.h>
#include <arpa/inet.h>
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
  if (stream->data_size == 0)
    stream->data_size = probe_data_size(data);
  ASSERT(stream->data_size == probe_data_size(data));
  size_t data_name_len = strlen(stream->data_name);
  size_t buff_size = sizeof(char)*header_len + 
                     sizeof(int) + sizeof(int) + 
                     sizeof(char)*data_name_len + 
                     sizeof(real_t) * (1 + stream->data_size);

  // Assemble the buffer.
  char buff[buff_size];

  // 1. The string "polymec-probe-stream".
  memcpy(buff, header, sizeof(char) * header_len); 
  size_t offset = header_len;

  // 2. An int containing the length of the probe's data name.
  int name_len = (int)data_name_len;
  memcpy(&buff[offset], &name_len, sizeof(int));
  offset += sizeof(int);

  // 3. An int containing the number of real numbers in the data.
  memcpy(&buff[offset], &(stream->data_size), sizeof(int));
  offset += sizeof(int);

  // 4. The characters in the probe's data name.
  memcpy(&buff[offset], stream->data_name, sizeof(char) * data_name_len);
  offset += sizeof(char) * data_name_len;

  // 5. A real_t containing the time of the acquisition.
  memcpy(&buff[offset], &t, sizeof(real_t));
  offset += sizeof(real_t);

  // 6. A sequence of real_t data.
  memcpy(&buff[offset], data->data, sizeof(real_t) * stream->data_size);
  offset += sizeof(real_t) * stream->data_size;

  // Write the buffer to the stream.
  write(stream->socket_fd, buff, buff_size);
}

bool probe_stream_on_acquire(probe_t* probe, const char* destination, int port)
{
  stream_context_t* context = polymec_malloc(sizeof(stream_context_t));
  context->data_name = probe_data_name(probe);
  context->data_size = 0;
  const char* p_name = probe_name(probe);

  // If we're given a valid port, try connecting via UDP first.
  char err_msg[129];
  if (port > 0)
  {
    context->socket_fd = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (context->socket_fd > 0)
    {
      memset(&(context->inet_addr), 0, sizeof(context->inet_addr));
      context->inet_addr.sin_family = AF_INET;
      inet_pton(AF_INET, destination, &context->inet_addr.sin_addr);
//      memcpy(&(context->inet_addr.sin_addr), destination, sizeof(char) * strlen(destination));
      context->inet_addr.sin_port = htons(port);
      int result = connect(context->socket_fd, (struct sockaddr*)&(context->inet_addr), sizeof(context->inet_addr));
      if (result != -1)
      {
        log_info("Probe %s: streaming %s to %s:%d via UDP", p_name, context->data_name, destination, port);
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
      log_info("Probe %s: streaming %s to %s via UNIX sockets", p_name, context->data_name, destination);
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
  log_info("Probe %s: Can't stream %s to %s (%s)", p_name, context->data_name, destination, err_msg);
  polymec_free(context);
  return false;
}

