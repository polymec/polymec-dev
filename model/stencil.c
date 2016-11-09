// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
                       int num_ghosts, exchanger_t* ex)
{
  ASSERT(num_indices > 0);
  ASSERT(offsets != NULL);
  ASSERT(indices != NULL);
  ASSERT(num_ghosts >= 0);
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  s->name = string_dup(name);
  s->num_indices = num_indices;
  s->offsets = offsets;
  s->indices = indices;
  s->weights = weights;
  s->num_ghosts = num_ghosts;
  s->ex = ex;
  return s;
}

stencil_t* unweighted_stencil_new(const char* name, int num_indices, 
                                  int* offsets, int* indices, 
                                  int num_ghosts, exchanger_t* ex)
{
  return stencil_new(name, num_indices, offsets, indices, NULL, num_ghosts, ex);
}

void stencil_free(stencil_t* stencil)
{
  polymec_free(stencil->name);
  polymec_free(stencil->offsets);
  polymec_free(stencil->indices);
  if (stencil->weights != NULL)
    polymec_free(stencil->weights);
  stencil->ex = NULL;
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
                     indices, weights, stencil->num_ghosts, ex);
}

void stencil_set_weights(stencil_t* stencil, real_t* weights)
{
  ASSERT(weights != NULL);
  if (stencil->weights != NULL)
    polymec_free(stencil->weights);
  stencil->weights = weights;
}

#if 0
stencil_t* augmented_stencil(stencil_t* stencil)
{
  MPI_Comm comm = exchanger_comm(stencil->ex);
  int rank, nproc;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nproc);

  // Gather the local neighbors of neighbors. 
  // n_of_n[i]->data[3*j] == rank owning the jth neighbor of i.
  // n_of_n[i]->data[3*j+1] == index of jth neighbor of i.
  // n_of_n[i]->data[3*j+2] == weight associated with jth neighbor of i.
  int N_local = stencil->num_indices;
  real_array_t* n_of_n[N_local]; 
  int_unordered_set_t* orig_neighbors = int_unordered_set_new();
  int max_index = N_local - 1;
  int max_nn = 0;
  for (int i = 0; i < N_local; ++i)
  {
    n_of_n[i] = real_array_new();

    int posj = 0, j;
    real_t wj;
    while (stencil_next(stencil, i, &posj, &j, &wj))
    {
      real_array_append(n_of_n[i], 1.0*rank);
      max_index = MAX(max_index, j);
      real_array_append(n_of_n[i], 1.0*j);
      real_array_append(n_of_n[i], wj);
      int_unordered_set_insert(orig_neighbors, j);
      int posk = 0, k;
      real_t wk;
      while (stencil_next(stencil, i, &posk, &k, &wk))
      {
        if (!int_unordered_set_contains(orig_neighbors, k))
        {
          max_index = MAX(max_index, k);
          real_array_append(n_of_n[i], rank);
          real_array_append(n_of_n[i], k);
          real_array_append(n_of_n[i], wk);
          int_unordered_set_insert(orig_neighbors, k);
        }
      }
    }

    max_nn = MAX(max_nn, n_of_n[i]->size);
    int_unordered_set_clear(orig_neighbors);
  }

  // Exchange n_of_n.
  serializer_t* s = int_array_serializer();
  exchanger_serialized_exchange(stencil->ex, n_of_n, 0, s);

  // At this point, we can read off all of the neighbors of neighbors of 
  // each of the indices, as well as the processes that own them.

  // Create an offsets array by inspecting n_of_n, and create the 
  // exchanger for the augmented stencil.
  int* offsets = polymec_malloc(sizeof(int) * (N_local+1));
  offsets[0] = N_local;
  exchanger_t* ex = exchanger_new(comm);
  int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* recv_map = int_ptr_unordered_map_new();
  for (int i = 0; i < N_local; ++i)
  {
    for (int j = 0; j < n_of_n[i]->size; ++j)
    {
      int proc = (int)n_of_n[i]->data[3*j];
      if (proc != rank)
      {
        int k = (int)n_of_n[i]->data[3*j+1];
        int* recv_p = int_ptr_unordered_map_get(recv_map, proc);
        int_array_t* recv;
        if (recv_p == NULL)
        {
          recv = int_array_new();
          int_ptr_unordered_map_insert_with_v_dtor(recv_map, proc, recv);
        }
        else
          recv = *recv_p;
        int_array_append(recv, k);
        ++(offsets[i+1]);
      }
    }
  }

  // Construct the augmented stencil and its exchanger.
  int* indices = polymec_malloc(sizeof(int) * offsets[nproc]);
  for (int i = 0; i < N_local; ++i)
  {
    for (int j = 0; j < n_of_n[i]->size; ++j)
    {
      int proc = (int)n_of_n[i]->data[3*j];
      if (proc != rank)
      {
        int k = (int)n_of_n[i]->data[3*j+1];
        real_t w = n_of_n[i]->data[3*j+2];
      }
    }
  }

  // Handle weights if needed.
  real_t* weights = NULL;
  if (stencil->weights != NULL)
  {
  }

  // Clean up.
  for (int i = 0; i < stencil->num_indices; ++i)
    real_array_free(n_of_n[i]);
  int_unordered_set_free(orig_neighbors);
}
#endif

static inline void swap_indices(int* i, int* j)
{
  int k = *i;
  *i = *j;
  *j = k;
}

static inline void swap_weights(real_t* wi, real_t* wj)
{
  real_t wk = *wi;
  *wi = *wj;
  *wj = wk;
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
          if (stencil->weights != NULL)
            swap_weights(&stencil->weights[j], &stencil->weights[last_index]);
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
  int num_indices, max_offset, num_weights, num_ghosts;
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
  byte_array_read_ints(bytes, 1, &num_ghosts, offset);

  // Exchanger stuff.
  serializer_t* ser = exchanger_serializer();
  exchanger_t* ex = serializer_read(ser, bytes, offset);

  return stencil_new(name, num_indices, offsets, indices, weights, num_ghosts, ex);
}

static void stencil_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  stencil_t* stencil = obj;

  // Write the name.
  int name_len = (int)strlen(stencil->name);
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
  byte_array_write_ints(bytes, 1, &(stencil->num_ghosts), offset);

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
  int N_local = stencil_num_indices(stencil), num_verts[nproc];
#if POLYMEC_HAVE_MPI
  MPI_Allgather(&N_local, 1, MPI_INT, num_verts, 1, MPI_INT, comm);
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

matrix_sparsity_t* sparsity_from_stencil(stencil_t* stencil)
{
  adj_graph_t* g = stencil_as_graph(stencil);
  return matrix_sparsity_from_graph(g, stencil->ex);
}

void silo_file_write_stencil(silo_file_t* file,
                             const char* stencil_name,
                             stencil_t* stencil)
{
  char name_name[FILENAME_MAX];
  snprintf(name_name, FILENAME_MAX, "%s_stencil_name", stencil_name);
  silo_file_write_string(file, name_name, stencil->name);

  char offsets_name[FILENAME_MAX];
  snprintf(offsets_name, FILENAME_MAX, "%s_stencil_offsets", stencil_name);
  silo_file_write_int_array(file, offsets_name, stencil->offsets, stencil->num_indices+1);

  char indices_name[FILENAME_MAX];
  snprintf(indices_name, FILENAME_MAX, "%s_stencil_indices", stencil_name);
  silo_file_write_int_array(file, indices_name, stencil->indices, stencil->offsets[stencil->num_indices]);

  char weights_name[FILENAME_MAX];
  snprintf(weights_name, FILENAME_MAX, "%s_stencil_weights", stencil_name);
  if (stencil->weights != NULL)
    silo_file_write_real_array(file, weights_name, stencil->weights, stencil->offsets[stencil->num_indices]);
  else
    silo_file_write_real_array(file, weights_name, NULL, 0);

  if (stencil->ex != NULL)
  {
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_stencil_ex", stencil_name);
    silo_file_write_exchanger(file, weights_name, stencil->ex);
  }
}

stencil_t* silo_file_read_stencil(silo_file_t* file,
                                  const char* stencil_name,
                                  MPI_Comm comm)
{
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  char name_name[FILENAME_MAX];
  snprintf(name_name, FILENAME_MAX, "%s_stencil_name", stencil_name);
  s->name = silo_file_read_string(file, name_name);

  char offsets_name[FILENAME_MAX];
  snprintf(offsets_name, FILENAME_MAX, "%s_stencil_offsets", stencil_name);
  size_t size;
  s->offsets = silo_file_read_int_array(file, offsets_name, &size);
  s->num_indices = (int)(size) - 1;

  char indices_name[FILENAME_MAX];
  snprintf(indices_name, FILENAME_MAX, "%s_stencil_indices", stencil_name);
  s->indices = silo_file_read_int_array(file, indices_name, &size);
  ASSERT((int)size == s->offsets[s->num_indices]);

  char weights_name[FILENAME_MAX];
  snprintf(weights_name, FILENAME_MAX, "%s_stencil_weights", stencil_name);
  size_t num_weights;
  s->weights = silo_file_read_real_array(file, weights_name, &num_weights);
  ASSERT((num_weights == s->offsets[s->num_indices]) || 
         ((num_weights == 0) && (s->weights == NULL)));

  char ex_name[FILENAME_MAX];
  snprintf(ex_name, FILENAME_MAX, "%s_stencil_ex", stencil_name);
  s->ex = silo_file_read_exchanger(file, ex_name, comm);
  return s;
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

  // Create an unweighted stencil.
  stencil_t* stencil = 
    unweighted_stencil_new("Distance-based point stencil", points->num_points, 
                           offsets, indices->data, num_ghosts, ex);

  // Clean up.
  int_array_release_data_and_free(indices); // Release control of index data.
  kd_tree_free(tree);

  // Add the ghost points to the point cloud if needed.
  point_cloud_set_num_ghosts(points, MAX(points->num_ghosts, num_ghosts));

  return stencil;
}
