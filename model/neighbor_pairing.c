// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "model/neighbor_pairing.h"

neighbor_pairing_t* neighbor_pairing_new(const char* name, int num_pairs, 
                                         int* pairs, real_t* weights,
                                         exchanger_t* ex)
{
  ASSERT(pairs != NULL);
  neighbor_pairing_t* p = polymec_malloc(sizeof(neighbor_pairing_t));
  p->name = string_dup(name);
  p->num_pairs = num_pairs;
  p->pairs = pairs;
  p->weights = weights;
  p->ex = ex;
  return p;
}

neighbor_pairing_t* unweighted_neighbor_pairing_new(const char* name, int num_pairs, 
                                                    int* pairs, exchanger_t* ex)
{
  return neighbor_pairing_new(name, num_pairs, pairs, NULL, ex);
}

void neighbor_pairing_free(neighbor_pairing_t* pairing)
{
  polymec_free(pairing->name);
  polymec_free(pairing->pairs);
  if (pairing->weights != NULL)
    polymec_free(pairing->weights);
  exchanger_free(pairing->ex);
  polymec_free(pairing);
}

void neighbor_pairing_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type)
{
  exchanger_exchange(pairing->ex, data, stride, tag, type);
}

int neighbor_pairing_start_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type)
{
  return exchanger_start_exchange(pairing->ex, data, stride, tag, type);
}

void neighbor_pairing_finish_exchange(neighbor_pairing_t* pairing, int token)
{
  exchanger_finish_exchange(pairing->ex, token);
}

static size_t np_byte_size(void* obj)
{
  neighbor_pairing_t* np = obj;
  
  // Data.
  size_t basic_storage = sizeof(int) + sizeof(char) * strlen(np->name) + 
                         sizeof(int) * (1 + 2*np->num_pairs + 1);
  if (np->weights != NULL)
    basic_storage += sizeof(real_t) * np->num_pairs;
  
  // Exchanger-related storage.
  serializer_t* ex_s = exchanger_serializer();
  size_t ex_storage = serializer_size(ex_s, np->ex);
  ex_s = NULL;

  return basic_storage + ex_storage;
}

static void* np_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the name.
  int name_len; 
  byte_array_read_ints(bytes, 1, &name_len, offset);
  char name[name_len+1];
  byte_array_read_chars(bytes, name_len, name, offset);
  name[name_len] = '\0';

  // Read the offsets, indices, weights.
  int num_pairs, num_weights;
  byte_array_read_ints(bytes, 1, &num_pairs, offset);
  int* pairs = polymec_malloc(sizeof(int) * 2 * num_pairs);
  byte_array_read_ints(bytes, 2*num_pairs, pairs, offset);
  byte_array_read_ints(bytes, 1, &num_weights, offset);
  real_t* weights = NULL;
  if (num_weights > 0)
    byte_array_read_reals(bytes, num_weights, weights, offset);

  // Exchanger stuff.
  serializer_t* ser = exchanger_serializer();
  exchanger_t* ex = serializer_read(ser, bytes, offset);

  return neighbor_pairing_new(name, num_pairs, pairs, weights, ex);
}

static void np_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  neighbor_pairing_t* np = obj;

  // Write the name.
  int name_len = strlen(np->name);
  byte_array_write_ints(bytes, 1, &name_len, offset);
  byte_array_write_chars(bytes, name_len, np->name, offset);

  // Write the offsets, indices, weights.
  byte_array_write_ints(bytes, 1, &np->num_pairs, offset);
  byte_array_write_ints(bytes, 2*np->num_pairs, np->pairs, offset);
  if (np->weights != NULL)
  {
    byte_array_write_ints(bytes, 1, &np->num_pairs, offset);
    byte_array_write_reals(bytes, np->num_pairs, np->weights, offset);
  }
  else
  {
    int zero = 0;
    byte_array_write_ints(bytes, 1, &zero, offset);
  }

  // Exchanger.
  serializer_t* ser = exchanger_serializer();
  serializer_write(ser, np->ex, bytes, offset);
  ser = NULL;
}

serializer_t* neighbor_pairing_serializer()
{
  return serializer_new(np_byte_size, np_byte_read, np_byte_write);
}

adj_graph_t* graph_from_point_cloud_and_neighbors(point_cloud_t* points, 
                                                  neighbor_pairing_t* neighbors)
{
  // Create a graph whose vertices are the cloud's points.
  int rank, nproc;
  MPI_Comm_size(points->comm, &nproc);
  MPI_Comm_rank(points->comm, &rank);
  adj_graph_t* g = adj_graph_new(points->comm, points->num_points);

  // Allocate space in the graph for the edges (neighbors associating points).
  int num_points = points->num_points;
  int* num_edges = polymec_malloc(sizeof(int) * num_points);
  memset(num_edges, 0, sizeof(int) * num_points);
  int pos = 0, i, j;
  while (neighbor_pairing_next(neighbors, &pos, &i, &j, NULL))
  {
    if (i < num_points)
      ++num_edges[i];
    if (j < num_points)
      ++num_edges[j];
  }
  for (int i = 0; i < num_points; ++i)
    adj_graph_set_num_edges(g, i, num_edges[i]);

  // Now fill in the edges.
  memset(num_edges, 0, sizeof(int) * num_points);
  pos = 0;
  while (neighbor_pairing_next(neighbors, &pos, &i, &j, NULL))
  {
    if (i < num_points)
    {
      int* edges = adj_graph_edges(g, i);
      edges[num_edges[i]++] = j;
    }
    if (j < num_points)
    {
      int* edges = adj_graph_edges(g, j);
      edges[num_edges[j]++] = i;
    }
  }

  // Clean up.
  polymec_free(num_edges);

  return g;
}

void silo_file_write_neighbor_pairing(silo_file_t* file,
                                      const char* neighbors_name,
                                      neighbor_pairing_t* neighbors)
{
  char name_name[FILENAME_MAX];
  snprintf(name_name, FILENAME_MAX, "%s_neighbor_pairing_name", neighbors_name);
  silo_file_write_string(file, name_name, neighbors->name);
  char pairs_name[FILENAME_MAX];
  snprintf(pairs_name, FILENAME_MAX, "%s_neighbor_pairing_pairs", neighbors_name);
  silo_file_write_int_array(file, pairs_name, neighbors->pairs, 2*neighbors->num_pairs);
  char weights_name[FILENAME_MAX];
  snprintf(weights_name, FILENAME_MAX, "%s_neighbor_pairing_weights", neighbors_name);
  if (neighbors->weights != NULL)
    silo_file_write_real_array(file, weights_name, neighbors->weights, neighbors->num_pairs);
  else
    silo_file_write_real_array(file, weights_name, NULL, 0);

  if (neighbors->ex != NULL)
  {
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_neighbor_pairing_ex", neighbors_name);
    silo_file_write_exchanger(file, weights_name, neighbors->ex);
  }
}

neighbor_pairing_t* silo_file_read_neighbor_pairing(silo_file_t* file,
                                                    const char* neighbors_name,
                                                    MPI_Comm comm)
{
  neighbor_pairing_t* p = polymec_malloc(sizeof(neighbor_pairing_t));
  char name_name[FILENAME_MAX];
  snprintf(name_name, FILENAME_MAX, "%s_neighbor_pairing_name", neighbors_name);
  p->name = silo_file_read_string(file, name_name);
  char pairs_name[FILENAME_MAX];
  snprintf(pairs_name, FILENAME_MAX, "%s_neighbor_pairing_pairs", neighbors_name);
  int size;
  p->pairs = silo_file_read_int_array(file, pairs_name, &size);
  ASSERT((size % 2) == 0);
  p->num_pairs = size/2;
  char weights_name[FILENAME_MAX];
  snprintf(weights_name, FILENAME_MAX, "%s_neighbor_pairing_weights", neighbors_name);
  int num_weights;
  p->weights = silo_file_read_real_array(file, weights_name, &num_weights);
  ASSERT((num_weights == p->num_pairs) || 
         ((num_weights == 0) && (p->weights == NULL)));
  char ex_name[FILENAME_MAX];
  snprintf(ex_name, FILENAME_MAX, "%s_neighbor_pairing_ex", neighbors_name);
  p->ex = silo_file_read_exchanger(file, ex_name, comm);
  return p;
}

