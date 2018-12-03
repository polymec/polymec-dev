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
  stream_context_t* stream = context;
  if (stream->data_size == 0)
    stream->data_size = probe_data_size(data);
  ASSERT(stream->data_size == probe_data_size(data));
  ASSERT(stream->data_size > 0);

  // Construct the string representation of the data.
  char data_str[20 * stream->data_size];
  data_str[0] = '[';
  size_t offset = 1;
  for (size_t i = 0; i < stream->data_size-1; ++i)
  {
    int len = sprintf(&data_str[offset], "%g,", data->data[i]);
    offset += len;
  }
  sprintf(&data_str[offset], "%g]%c", data->data[stream->data_size-1], '\0');

  // Figure out the size of the JSON payload.
  const char* json_format_str = "{\"name\": \"%s\", \"time\": %g, \"data\": %s}";
  size_t json_size = sizeof(char) * 
                     (strlen(json_format_str) + strlen(stream->data_name) + \
                      20 + strlen(data_str));

  // Assemble the JSON object text representation.
  char json[json_size];
  sprintf(json, json_format_str, stream->data_name, t, data_str);

  // Write the buffer to the stream.
  write(stream->socket_fd, json, strlen(json));
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
  else
    snprintf(err_msg, 128, "Invalid port (%d).", port);

  // Better luck next time!
  log_info("Probe %s: Can't stream %s to %s (%s)", p_name, context->data_name, destination, err_msg);
  polymec_free(context);
  return false;
}

