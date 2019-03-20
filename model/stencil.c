// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/stencil.h"
#include "core/timer.h"
#include "core/unordered_set.h"
#include "core/array.h"
#include "core/kd_tree.h"

stencil_t* stencil_new(const char* name, size_t num_indices,
                       int* offsets, int* indices,
                       size_t num_ghosts, exchanger_t* ex)
{
  ASSERT(num_indices > 0);
  ASSERT(offsets != NULL);
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  s->name = string_dup(name);
  s->num_indices = num_indices;
  s->offsets = offsets;
  s->indices = indices;
  s->num_ghosts = num_ghosts;
  s->ex = ex;
  return s;
}

void stencil_free(stencil_t* stencil)
{
  polymec_free(stencil->name);
  polymec_free(stencil->offsets);
  if (stencil->indices != NULL)
    polymec_free(stencil->indices);
  stencil->ex = NULL;
  polymec_free(stencil);
}

stencil_t* stencil_clone(stencil_t* stencil)
{
  int* offsets = polymec_malloc(sizeof(int) * (stencil->num_indices+1));
  memcpy(offsets, stencil->offsets, sizeof(int) * (stencil->num_indices+1));
  int size = offsets[stencil->num_indices];
  int* indices = polymec_malloc(sizeof(int) * size);
  memcpy(indices, stencil->indices, sizeof(int) * size);
  exchanger_t* ex = exchanger_clone(stencil->ex);
  return stencil_new(stencil->name, stencil->num_indices, offsets,
                     indices, stencil->num_ghosts, ex);
}

static inline void swap_indices(int* i, int* j)
{
  int k = *i;
  *i = *j;
  *j = k;
}

void stencil_trim(stencil_t* stencil, int_unordered_set_t** neighbors_to_trim)
{
  int num_trimmed_indices[stencil->num_indices];
  memset(num_trimmed_indices, 0, sizeof(int) * stencil->num_indices);

  for (int i = 0; i < stencil->num_indices; ++i)
  {
    int_unordered_set_t* indices = neighbors_to_trim[i];
    if ((indices != NULL) && (indices->size > 0))
    {
      int start = stencil->offsets[i];
      int end = stencil->offsets[i+1];
      for (int j = start; j < end; ++j)
      {
        int k = stencil->indices[j];
        if (int_unordered_set_contains(indices, k))
        {
          // Remove k from i's neighbor list by swapping it to the back.
          int last_index = end - 1 - num_trimmed_indices[i];
          swap_indices(&stencil->indices[j], &stencil->indices[last_index]);
          ++num_trimmed_indices[i];
          if (num_trimmed_indices[i] == indices->size)
            break;
        }
      }
    }
  }

  // Now re-index the neighbors.
  for (int i = 0; i < stencil->num_indices; ++i)
  {
    int decrease = 0;
    for (int j = 0; j < i; ++j)
      decrease += num_trimmed_indices[j];
    if (decrease > 0)
    {
      int src = stencil->offsets[i];
      int dest = stencil->offsets[i] - decrease;
      int num_indices = stencil->offsets[i+1] - stencil->offsets[i] - num_trimmed_indices[i];
      memmove(&stencil->indices[dest], &stencil->indices[src],
              sizeof(int) * num_indices);
      stencil->offsets[i] -= decrease;
    }
  }
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
                         sizeof(int) * (1 + stencil->num_indices + 1 + stencil->offsets[stencil->num_indices] + 1 + 1);

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

  // Read the offsets, indices.
  size_t num_indices, num_ghosts;
  int max_offset;
  byte_array_read_size_ts(bytes, 1, &num_indices, offset);
  int* offsets = polymec_malloc(sizeof(int) * (num_indices+1));
  byte_array_read_ints(bytes, num_indices+1, offsets, offset);
  byte_array_read_ints(bytes, 1, &max_offset, offset);
  int* indices = polymec_malloc(sizeof(int) * max_offset);
  byte_array_read_ints(bytes, max_offset, indices, offset);
  byte_array_read_size_ts(bytes, 1, &num_ghosts, offset);

  // Exchanger stuff.
  serializer_t* ser = exchanger_serializer();
  exchanger_t* ex = serializer_read(ser, bytes, offset);

  return stencil_new(name, num_indices, offsets, indices, num_ghosts, ex);
}

static void stencil_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  stencil_t* stencil = obj;

  // Write the name.
  int name_len = (int)strlen(stencil->name);
  byte_array_write_ints(bytes, 1, &name_len, offset);
  byte_array_write_chars(bytes, name_len, stencil->name, offset);

  // Write the offsets, indices.
  byte_array_write_size_ts(bytes, 1, &stencil->num_indices, offset);
  byte_array_write_ints(bytes, stencil->num_indices+1, stencil->offsets, offset);
  int max_offset = stencil->offsets[stencil->num_indices];
  byte_array_write_ints(bytes, 1, &max_offset, offset);
  byte_array_write_ints(bytes, max_offset, stencil->indices, offset);
  byte_array_write_size_ts(bytes, 1, &(stencil->num_ghosts), offset);

  // Exchanger.
  serializer_t* ser = exchanger_serializer();
  serializer_write(ser, stencil->ex, bytes, offset);
  ser = NULL;
}

adj_graph_t* stencil_as_graph(stencil_t* stencil)
{
  // Create a graph whose vertices are the cloud's points.
  MPI_Comm comm = exchanger_comm(stencil->ex);
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  size_t N_local = stencil_num_indices(stencil), num_verts[nproc];
#if POLYMEC_HAVE_MPI
  MPI_Allgather(&N_local, 1, MPI_SIZE_T, num_verts, 1, MPI_SIZE_T, comm);
#else
  num_verts[0] = N_local;
#endif
  index_t vtx_dist[nproc+1];
  vtx_dist[0] = 0;
  for (int p = 0; p < nproc; ++p)
    vtx_dist[p+1] = vtx_dist[p] + num_verts[p];
  return adj_graph_from_arrays(comm, vtx_dist,
                               stencil->indices, stencil->offsets,
                               false);
}

serializer_t* stencil_serializer()
{
  return serializer_new("stencil", stencil_byte_size, stencil_byte_read, stencil_byte_write, DTOR(stencil_free));
}

matrix_sparsity_t* matrix_sparsity_from_stencil(stencil_t* stencil)
{
  adj_graph_t* g = stencil_as_graph(stencil);
  matrix_sparsity_t* sp = matrix_sparsity_from_graph(g, stencil->ex);
  adj_graph_free(g);
  return sp;
}

stencil_t* distance_based_point_stencil_new(point_cloud_t* points,
                                            real_t* R)
{
  ASSERT(R != NULL);
#ifndef NDEBUG
  for (int i = 0; i < points->num_points; ++i)
    ASSERT(R[i] > 0.0);
#endif

  // Stick all the points into a kd-tree so that we can pair them up.
  kd_tree_t* tree = kd_tree_new(points->points, points->num_points);

  // Find the maximum radius of interaction.
  real_t R_max = -REAL_MAX;
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
    // Find all the neighbors for this point.
    point_t* xi = &points->points[i];
    int_array_t* neighbors = kd_tree_within_radius(tree, xi, R[i]);
    offsets[i+1] = offsets[i] + (int)neighbors->size;
    for (int k = 0; k < neighbors->size; ++k)
    {
      int j = neighbors->data[k];
      int_array_append(indices, j);
    }
    int_array_free(neighbors);
  }

  // Find the number of ghost points referred to within the stencil.
  int num_ghosts = (int)(kd_tree_size(tree) - points->num_points);

  // Create the stencil.
  stencil_t* stencil =
    stencil_new("Distance-based point stencil", points->num_points,
                offsets, indices->data, num_ghosts, ex);

  // Clean up.
  int_array_release_data_and_free(indices); // Release control of index data.
  kd_tree_free(tree);

  // Add the ghost points to the point cloud if needed.
  point_cloud_set_num_ghosts(points, MAX(points->num_ghosts, num_ghosts));

  return stencil;
}
