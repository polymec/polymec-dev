// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/stencil.h"
#include "core/unordered_set.h"
#include "core/array.h"
#include "core/kd_tree.h"

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

stencil_t* stencil_clone(stencil_t* stencil)
{
  int* offsets = polymec_malloc(sizeof(int) * (stencil->num_indices+1));
  memcpy(offsets, stencil->offsets, sizeof(real_t) * (stencil->num_indices+1));
  int size = offsets[stencil->num_indices];
  int* indices = polymec_malloc(sizeof(int) * size);
  memcpy(indices, stencil->indices, sizeof(real_t) * size);
  real_t* weights = NULL;
  if (stencil->weights != NULL)
  {
    weights = polymec_malloc(sizeof(real_t) * size);
    memcpy(weights, stencil->weights, sizeof(real_t) * size);
  }
  exchanger_t* ex = exchanger_clone(stencil->ex);
  return stencil_new(string_dup(stencil->name), stencil->num_indices, offsets, 
                     indices, weights, ex);
}

void stencil_set_weights(stencil_t* stencil, real_t* weights)
{
  ASSERT(weights != NULL);
  if (stencil->weights != NULL)
    polymec_free(stencil->weights);
  stencil->weights = weights;
}

void stencil_augment(stencil_t* stencil)
{
  // Gather the local neighbors of neighbors.
  int_array_t* n_of_n[stencil->num_indices]; 
  real_array_t* n_of_n_w[stencil->num_indices]; 
  int_unordered_set_t* old_neighbors = int_unordered_set_new();
  int max_index = stencil->num_indices - 1;
  int max_nn = 0;
  for (int i = 0; i < stencil->num_indices; ++i)
  {
    n_of_n[i] = int_array_new();
    n_of_n_w[i] = real_array_new();

    int posj = 0, j;
    real_t wj;
    while (stencil_next(stencil, i, &posj, &j, &wj))
    {
      max_index = MAX(max_index, j);
      int_array_append(n_of_n[i], j);
      real_array_append(n_of_n_w[i], wj);
      int_unordered_set_insert(old_neighbors, j);
      int posk = 0, k;
      real_t wk;
      while (stencil_next(stencil, i, &posk, &k, &wk))
      {
        if (!int_unordered_set_contains(old_neighbors, k))
        {
          max_index = MAX(max_index, k);
          int_array_append(n_of_n[i], k);
          real_array_append(n_of_n_w[i], wk);
          int_unordered_set_insert(old_neighbors, k);
        }
      }
    }

    max_nn = MAX(max_nn, n_of_n[i]->size);
    int_unordered_set_clear(old_neighbors);
  }

  // Exchange the number of neighbors for interior indices.
  int nn[stencil->num_indices];
  for (int i = 0; i < stencil->num_indices; ++i)
    nn[i] = n_of_n[i]->size;
  exchanger_exchange(stencil->ex, nn, 1, 0, MPI_INT);

  // Set up a new exchanger.
  exchanger_t* ex = exchanger_new(exchanger_comm(stencil->ex));
  int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* recv_map = int_ptr_unordered_map_new();

  // Clean up.
  for (int i = 0; i < stencil->num_indices; ++i)
  {
    int_array_free(n_of_n[i]);
    real_array_free(n_of_n_w[i]);
  }
  int_unordered_set_free(old_neighbors);
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

stencil_t* distance_based_point_stencil_new(point_cloud_t* points,
                                            real_t* R,
                                            int* num_ghost_points)
{
  ASSERT(R != NULL);
  ASSERT(num_ghost_points != NULL);
#ifndef NDEBUG
  for (int i = 0; i < points->num_points; ++i)
    ASSERT(R[i] > 0.0);
#endif

  // Stick all the points into a kd-tree so that we can pair them up.
  kd_tree_t* tree = kd_tree_new(points->points, points->num_points);

  // Find the maximum radius of interaction.
  real_t R_max = -FLT_MAX;
  for (int i = 0; i < points->num_points; ++i)
    R_max = MAX(R_max, R[i]);

  // Add ghost points to the kd-tree and fetch an exchanger. This may add 
  // too many ghost points, but hopefully that won't be an issue.
  exchanger_t* ex = kd_tree_find_ghost_points(tree, points->comm, R_max);

  // We'll toss stencil data into these arrays.
  int* offsets = polymec_malloc(sizeof(int) * (points->num_points+1));
  int_array_t* indices = int_array_new();
  offsets[0] = 0;
  for (int i = 0; i < points->num_points; ++i)
  {
    // Find all the neighbors for this point. We only count those 
    // neighbors {j} for which j > i.
    point_t* xi = &points->points[i];
    int_array_t* neighbors = kd_tree_within_radius(tree, xi, R_max);
    offsets[i+1] = offsets[i] + neighbors->size;
    for (int k = 0; k < neighbors->size; ++k)
    {
      int j = neighbors->data[k];
      if (j > i)
      {
        real_t D = point_distance(xi, &points->points[j]);
        if (D < R[i])
          int_array_append(indices, j);
      }
    }
    int_array_free(neighbors);
  }

  // Create an unweighted stencil.
  stencil_t* stencil = 
    unweighted_stencil_new("Distance-based point stencil", points->num_points, 
                           offsets, indices->data, ex);

  // Set the number of ghost points referred to within the neighbor pairing.
  *num_ghost_points = kd_tree_size(tree) - points->num_points;

  // Clean up.
  int_array_release_data_and_free(indices); // Release control of index data.
  kd_tree_free(tree);

  return stencil;
}
