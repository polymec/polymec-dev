// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/blob_exchanger.h"

// A blob buffer is really just a big chunk of memory with an associated type.
struct blob_buffer_t 
{
  void* b;
};

void blob_buffer_free(blob_buffer_t* buffer)
{
  polymec_free(buffer->b);
  polymec_free(buffer);
}

struct blob_exchanger_t 
{
  MPI_Comm comm;
};

void blob_exchanger_proc_map_add_index(blob_exchanger_proc_map_t* map, 
                                       int process, 
                                       int blob_index)
{
}

static void blob_exchanger_free(void* context)
{
//  blob_exchanger_t* ex = context;
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

blob_buffer_t* blob_exchanger_create_buffer(blob_exchanger_t* ex)
{
  blob_buffer_t* b = polymec_malloc(sizeof(blob_buffer_t));
  return b;
}

void blob_exchanger_exchange(blob_exchanger_t* ex, 
                             int tag,
                             blob_buffer_t* buffer)
{
}

int blob_exchanger_start_exchange(blob_exchanger_t* ex, 
                                  int tag,
                                  blob_buffer_t* buffer)
{
  return 0; // FIXME
}

void blob_exchanger_finish_exchange(blob_exchanger_t* ex, int token)
{
}

bool blob_exchanger_next_send_blob(blob_exchanger_t* ex, 
                                   blob_buffer_t* buffer,
                                   int* pos, 
                                   int* remote_process, 
                                   int* blob_index, 
                                   size_t* blob_size)
{
  return false;
}

bool blob_exchanger_next_receive_blob(blob_exchanger_t* ex, 
                                      blob_buffer_t* buffer,
                                      int* pos, 
                                      int* remote_process, 
                                      int* blob_index, 
                                      size_t* blob_size)
{
  return false;
}

void blob_exchanger_copy_in(blob_exchanger_t* ex,
                            int blob_index,
                            void* blob,
                            blob_buffer_t* buffer)
{
}

void blob_exchanger_copy_out(blob_exchanger_t* ex,
                             blob_buffer_t* buffer,
                             int blob_index,
                             void* blob)
{
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
}

void blob_exchanger_disable_deadlock_detection(blob_exchanger_t* ex)
{
}

bool blob_exchanger_deadlock_detection_enabled(blob_exchanger_t* ex)
{
  return false;
}

void blob_exchanger_fprintf(blob_exchanger_t* ex, FILE* stream)
{
}

