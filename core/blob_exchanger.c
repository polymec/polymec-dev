// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/blob_exchanger.h"

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
  blob_exchanger_t* ex = context;
}

blob_exchanger_t* blob_exchanger_new(MPI_Comm comm,
                                     blob_exchanger_proc_map* send_map,
                                     blob_exchanger_proc_map* receive_map,
                                     blob_exchanger_size_map* blob_size_map)
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

void* blob_exchanger_new_send_buffer(blob_exchanger_t* ex)
{
}

void* blob_exchanger_new_receive_buffer(blob_exchanger_t* ex)
{
}

void blob_exchanger_exchange(blob_exchanger_t* ex, 
                             int tag,
                             void* send_buffer, 
                             void* receive_buffer)
{
}

int blob_exchanger_start_exchange(blob_exchanger_t* ex, 
                                  int tag,
                                  void* send_buffer, 
                                  void* receive_buffer)
{
}

void blob_exchanger_finish_exchange(blob_exchanger_t* ex, int token)
{
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
}

void blob_exchanger_fprintf(blob_exchanger_t* ex, FILE* stream)
{
}

bool blob_exchanger_next_send_blob(blob_exchanger_t* ex, 
                                   void* send_buffer,
                                   int* pos, 
                                   int* remote_process, 
                                   int* blob_index, 
                                   void** blob, 
                                   size_t* blob_size)
{
}

bool blob_exchanger_next_receive_blob(blob_exchanger_t* ex, 
                                      void* receive_buffer,
                                      int* pos, 
                                      int* remote_process, 
                                      int* blob_index, 
                                      void** blob, 
                                      size_t* blob_size)
{
}

bool blob_exchanger_verify(blob_exchanger_t* ex, 
                           void (*handler)(const char* format, ...))
{
}

