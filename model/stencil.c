// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/stencil.h"
#include "core/unordered_set.h"
#include "core/array.h"

stencil_t* stencil_new(const char* name, int num_indices, 
                       int* offsets, int* indices, real_t* weights,
                       exchanger_t* ex)
{
  ASSERT(num_indices > 0);
  ASSERT(offsets != NULL);
  ASSERT(indices != NULL);
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  s->name = string_dup(name);
  s->num_indices = num_indices;
  s->offsets = offsets;
  s->indices = indices;
  s->weights = weights;
  s->ex = ex;
  return s;
}

stencil_t* unweighted_stencil_new(const char* name, int num_indices, 
                                  int* offsets, int* indices, exchanger_t* ex)
{
  return stencil_new(name, num_indices, offsets, indices, NULL, ex);
}

void stencil_free(stencil_t* stencil)
{
  polymec_free(stencil->name);
  polymec_free(stencil->offsets);
  polymec_free(stencil->indices);
  if (stencil->weights != NULL)
    polymec_free(stencil->weights);
  exchanger_free(stencil->ex);
  polymec_free(stencil);
}

void stencil_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type)
{
  exchanger_exchange(stencil->ex, data, stride, tag, type);
}

int stencil_start_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type)
{
  return exchanger_start_exchange(stencil->ex, data, stride, tag, type);
}

void stencil_finish_exchange(stencil_t* stencil, int token)
{
  exchanger_finish_exchange(stencil->ex, token);
}

exchanger_t* stencil_exchanger(stencil_t* stencil)
{
  return stencil->ex;
}

static size_t stencil_byte_size(void* obj)
{
  stencil_t* stencil = obj;
  
  // Data.
  size_t basic_storage = sizeof(int) + sizeof(char) * strlen(stencil->name) + 
                         sizeof(int) * (1 + stencil->num_indices + 1 + stencil->offsets[stencil->num_indices] + 1);
  if (stencil->weights != NULL)
    basic_storage += sizeof(real_t) * stencil->offsets[stencil->num_indices];
  
  // Exchanger-related storage.
  serializer_t* ex_s = exchanger_serializer();
  size_t ex_storage = serializer_size(ex_s, stencil->ex);
  ex_s = NULL;

  return basic_storage + ex_storage;
}

static void* stencil_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the name.
  int name_len; 
  byte_array_read_ints(bytes, 1, &name_len, offset);
  char name[name_len+1];
  byte_array_read_chars(bytes, name_len, name, offset);
  name[name_len] = '\0';

  // Read the offsets, indices, weights.
  int num_indices, max_offset, num_weights;
  byte_array_read_ints(bytes, 1, &num_indices, offset);
  int* offsets = polymec_malloc(sizeof(int) * (num_indices+1));
  byte_array_read_ints(bytes, num_indices+1, offsets, offset);
  byte_array_read_ints(bytes, 1, &max_offset, offset);
  int* indices = polymec_malloc(sizeof(int) * max_offset);
  byte_array_read_ints(bytes, max_offset, indices, offset);
  byte_array_read_ints(bytes, 1, &num_weights, offset);
  real_t* weights = NULL;
  if (num_weights > 0)
    byte_array_read_real_ts(bytes, num_weights, weights, offset);

  // Exchanger stuff.
  serializer_t* ser = exchanger_serializer();
  exchanger_t* ex = serializer_read(ser, bytes, offset);

  return stencil_new(name, num_indices, offsets, indices, weights, ex);
}

static void stencil_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  stencil_t* stencil = obj;

  // Write the name.
  int name_len = strlen(stencil->name);
  byte_array_write_ints(bytes, 1, &name_len, offset);
  byte_array_write_chars(bytes, name_len, stencil->name, offset);

  // Write the offsets, indices, weights.
  byte_array_write_ints(bytes, 1, &stencil->num_indices, offset);
  byte_array_write_ints(bytes, stencil->num_indices+1, stencil->offsets, offset);
  int max_offset = stencil->offsets[stencil->num_indices];
  byte_array_write_ints(bytes, 1, &max_offset, offset);
  byte_array_write_ints(bytes, max_offset, stencil->indices, offset);
  if (stencil->weights != NULL)
  {
    byte_array_write_ints(bytes, 1, &max_offset, offset);
    byte_array_write_real_ts(bytes, max_offset, stencil->weights, offset);
  }
  else
  {
    int zero = 0;
    byte_array_write_ints(bytes, 1, &zero, offset);
  }

  // Exchanger.
  serializer_t* ser = exchanger_serializer();
  serializer_write(ser, stencil->ex, bytes, offset);
  ser = NULL;
}

serializer_t* stencil_serializer()
{
  return serializer_new("stencil", stencil_byte_size, stencil_byte_read, stencil_byte_write, DTOR(stencil_free));
}

adj_graph_t* graph_from_point_cloud_and_stencil(point_cloud_t* points, 
                                                stencil_t* stencil)
{
  ASSERT(stencil_num_indices(stencil) == points->num_points);

  // Create a graph whose vertices are the cloud's points.
  int rank, nproc;
  MPI_Comm_size(points->comm, &nproc);
  MPI_Comm_rank(points->comm, &rank);
  adj_graph_t* g = adj_graph_new(points->comm, points->num_points);

  // Allocate space in the graph for the edges (stencil size).
  int num_points = points->num_points;
  for (int i = 0; i < num_points; ++i)
    adj_graph_set_num_edges(g, i, stencil_size(stencil, i));

  // Now fill in the edges.
  for (int i = 0; i < num_points; ++i)
  {
    int* edges = adj_graph_edges(g, i);
    int offset = 0, pos = 0, j;
    while (stencil_next(stencil, i, &pos, &j, NULL))
      edges[offset++] = j;
  }

  return g;
}

