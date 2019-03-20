// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/kd_tree.h"
#include "core/array.h"
#include "core/unordered_set.h"
#include "model/neighbor_pairing.h"

neighbor_pairing_t* neighbor_pairing_new(const char* name, size_t num_pairs,
                                         int* pairs, exchanger_t* ex)
{
  ASSERT(pairs != NULL);
  ASSERT(ex != NULL);
  neighbor_pairing_t* p = polymec_malloc(sizeof(neighbor_pairing_t));
  p->name = string_dup(name);
  p->num_pairs = num_pairs;
  p->pairs = pairs;
  p->ex = ex;
  return p;
}

void neighbor_pairing_free(neighbor_pairing_t* pairing)
{
  polymec_free(pairing->name);
  polymec_free(pairing->pairs);
  pairing->ex = NULL;
  polymec_free(pairing);
}

exchanger_t* neighbor_pairing_exchanger(neighbor_pairing_t* pairing)
{
  return pairing->ex;
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

  // Read the offsets, indices.
  size_t num_pairs;
  byte_array_read_size_ts(bytes, 1, &num_pairs, offset);
  int* pairs = polymec_malloc(sizeof(int) * 2 * num_pairs);
  byte_array_read_ints(bytes, 2*num_pairs, pairs, offset);

  // Exchanger stuff.
  serializer_t* ser = exchanger_serializer();
  exchanger_t* ex = serializer_read(ser, bytes, offset);

  return neighbor_pairing_new(name, num_pairs, pairs, ex);
}

static void np_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  neighbor_pairing_t* np = obj;

  // Write the name.
  int name_len = (int)strlen(np->name);
  byte_array_write_ints(bytes, 1, &name_len, offset);
  byte_array_write_chars(bytes, name_len, np->name, offset);

  // Write the offsets, indices.
  byte_array_write_size_ts(bytes, 1, &np->num_pairs, offset);
  byte_array_write_ints(bytes, 2*np->num_pairs, np->pairs, offset);

  // Exchanger.
  serializer_t* ser = exchanger_serializer();
  serializer_write(ser, np->ex, bytes, offset);
  ser = NULL;
}

serializer_t* neighbor_pairing_serializer()
{
  return serializer_new("neighbor_pairing", np_byte_size, np_byte_read, np_byte_write, DTOR(neighbor_pairing_free));
}

static void free_pair(int* pair)
{
  polymec_free(pair);
}

neighbor_pairing_t* neighbor_pairing_from_stencil(stencil_t* stencil)
{
  int_array_t* pairs = int_array_new();
  int_pair_int_unordered_map_t* pair_map = int_pair_int_unordered_map_new();

  // Extract all the pairs from the stencil.
  size_t num_indices = stencil_num_indices(stencil);
  for (int i = 0; i < num_indices; ++i)
  {
    int pos = 0, j;
    while (stencil_next(stencil, i, &pos, &j))
    {
      int small = MIN(i, j);
      int big = MAX(i, j);
      int pair[2] = {small, big};
      int* pair_index_p = int_pair_int_unordered_map_get(pair_map, pair);
      if (pair_index_p == NULL)
      {
        int pair_index = (int)(pairs->size/2);
        int_array_append(pairs, small);
        int_array_append(pairs, big);
        int* p = polymec_malloc(sizeof(int) * 2);
        p[0] = small; p[1] = big;
        int_pair_int_unordered_map_insert_with_k_dtor(pair_map, p, pair_index, free_pair);
      }
    }
  }

  // The exchanger should be the same for both of these things.
  exchanger_t* ex = exchanger_clone(stencil->ex);

  // Build the neighbor pairing.
  neighbor_pairing_t* neighbors = neighbor_pairing_new(stencil->name,
                                                       pairs->size/2,
                                                       pairs->data,
                                                       ex);

  // Let the neighbor pairing steal the array data.
  int_array_release_data_and_free(pairs);
  int_pair_int_unordered_map_free(pair_map);

  return neighbors;
}

stencil_t* stencil_from_point_cloud_and_neighbors(point_cloud_t* points,
                                                  neighbor_pairing_t* neighbors)
{
  // Count up the numbers of neighbors for each index.
  int_int_unordered_map_t* counts = int_int_unordered_map_new();
  int pos = 0, i, j;
  while (neighbor_pairing_next(neighbors, &pos, &i, &j))
  {
    int small = MIN(i, j);
    int big = MAX(i, j);
    int* small_n_p = int_int_unordered_map_get(counts, small);
    if (small_n_p == NULL)
      int_int_unordered_map_insert(counts, small, 1);
    else
      ++(*small_n_p);
    int* big_n_p = int_int_unordered_map_get(counts, big);
    if (big_n_p == NULL)
      int_int_unordered_map_insert(counts, big, 1);
    else
      ++(*big_n_p);
  }

  // The exchanger should be the same for both of these things.
  exchanger_t* ex = exchanger_clone(neighbors->ex);

  // Set up the arrays for the stencil data.
  size_t num_indices = points->num_points;
  int* offsets = polymec_malloc(sizeof(int) * (num_indices+1));
  offsets[0] = 0;
  for (int k = 0; k < num_indices; ++k)
  {
    int* count_p = int_int_unordered_map_get(counts, k);
    int count = (count_p != NULL) ? *count_p : 0;
    offsets[k+1] = offsets[k] + count;
  }
  int N = offsets[num_indices];
  int* indices = polymec_malloc(sizeof(int) * N);
  int_int_unordered_map_free(counts);

  // Now extract the data from the neighbor pairing.
  size_t num_ghosts = points->num_ghosts;
  pos = 0;
  int which[num_indices + num_ghosts];
  memset(which, 0, sizeof(int) * (num_indices + num_ghosts));
  while (neighbor_pairing_next(neighbors, &pos, &i, &j))
  {
    if (i < num_indices)
    {
      indices[offsets[i]+which[i]] = j;
      ++which[i];
    }
    if (j < num_indices)
    {
      indices[offsets[j]+which[j]] = i;
      ++which[j];
    }
  }

  // Construct the stencil.
  return stencil_new(neighbors->name, num_indices,
                     offsets, indices, num_ghosts, ex);
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
  size_t num_points = points->num_points;
  int* num_edges = polymec_calloc(num_points, sizeof(int));
  {
    int pos = 0, i, j;
    while (neighbor_pairing_next(neighbors, &pos, &i, &j))
    {
      if (i < num_points)
        ++num_edges[i];
      if (j < num_points)
        ++num_edges[j];
    }
  }
  for (int i = 0; i < num_points; ++i)
    adj_graph_set_num_edges(g, i, num_edges[i]);

  // Now fill in the edges.
  memset(num_edges, 0, sizeof(int) * num_points);
  {
    int pos = 0, i, j;
    while (neighbor_pairing_next(neighbors, &pos, &i, &j))
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
  }

  // Clean up.
  polymec_free(num_edges);

  return g;
}

matrix_sparsity_t* sparsity_from_point_cloud_and_neighbors(point_cloud_t* points,
                                                           neighbor_pairing_t* neighbors)
{
  // Figure out the domain decomposition.
  MPI_Comm comm = points->comm;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  index_t num_points[nproc];
  index_t num_local_points = (index_t)points->num_points;
  MPI_Allgather(&num_local_points, 1, MPI_INDEX_T,
                num_points, 1, MPI_INDEX_T, comm);
  index_t row_dist[nproc+1];
  row_dist[0] = 0;
  for (int p = 0; p < nproc; ++p)
    row_dist[p+1] = row_dist[p] + num_points[p];

  // Get global indices for the locally-represented points.
  index_t global_ids[points->num_points + points->num_ghosts];
  for (int i = 0; i < points->num_points; ++i)
    global_ids[i] = (index_t)(row_dist[rank] + i);
  neighbor_pairing_exchange(neighbors, global_ids, 1, 0, MPI_INDEX_T);

  // Create a matrix sparsity pattern using the given neighbors and
  // allocate column space.
  matrix_sparsity_t* sparsity = matrix_sparsity_new(comm, row_dist);
  index_t num_cols[points->num_points];
  memset(num_cols, 0, sizeof(index_t) * points->num_points);
  int pos = 0, i, j;
  while (neighbor_pairing_next(neighbors, &pos, &i, &j))
  {
    ++num_cols[i];
    ++num_cols[j];
  }
  index_t offsets[points->num_points];
  for (int k = 0; k < points->num_points; ++k)
  {
    matrix_sparsity_set_num_columns(sparsity, global_ids[i], num_cols[i]);

    // Add the diagonal entry, while we're here.
    index_t* columns = matrix_sparsity_columns(sparsity, global_ids[k]);
    columns[0] = global_ids[k];
    offsets[k] = 1;
  }

  // Now step through and add each (i, j) pair.
  pos = 0;
  while (neighbor_pairing_next(neighbors, &pos, &i, &j))
  {
    index_t* i_columns = matrix_sparsity_columns(sparsity, global_ids[i]);
    index_t* j_columns = matrix_sparsity_columns(sparsity, global_ids[j]);
    i_columns[offsets[i]] = global_ids[j];
    ++offsets[i];
    j_columns[offsets[j]] = global_ids[i];
    ++offsets[j];
  }

  return sparsity;
}

neighbor_pairing_t* distance_based_neighbor_pairing_new(point_cloud_t* points,
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
  real_t R_max = -REAL_MAX;
  for (int i = 0; i < points->num_points; ++i)
    R_max = MAX(R_max, R[i]);

  // Add ghost points to the kd-tree and fetch an exchanger. This may add
  // too many ghost points, but hopefully that won't be an issue.
  exchanger_t* ex = kd_tree_find_ghost_points(tree, points->comm, R_max);

  // Make parallel-sensible coordinate and radius fields with ghost candidate
  // values filled in.
  size_t tree_size = kd_tree_size(tree);
  point_t* x_par = polymec_malloc(sizeof(point_t) * tree_size);
  memcpy(x_par, points->points, 3 * sizeof(real_t) * points->num_points);
  exchanger_exchange(ex, x_par, 3, 0, MPI_REAL_T);
  real_t* R_par = polymec_malloc(sizeof(real_t) * tree_size);
  memcpy(R_par, R, sizeof(real_t) * points->num_points);
  exchanger_exchange(ex, R_par, 1, 0, MPI_REAL_T);

  // We'll toss neighbor pairs into this expandable array.
  int_array_t* pair_array = int_array_new();

  for (int i = 0; i < points->num_points; ++i)
  {
    // Find all the neighbors for this point. We only count those
    // neighbors {j} for which j > i.
    point_t* xi = &x_par[i];
    int_array_t* neighbors = kd_tree_within_radius(tree, xi, R_max);
    for (int k = 0; k < neighbors->size; ++k)
    {
      int j = neighbors->data[k];
      if (j > i)
      {
        real_t D = point_distance(xi, &x_par[j]);
        if (D < MAX(R_par[i], R_par[j]))
        {
          int_array_append(pair_array, i);
          int_array_append(pair_array, j);
        }
      }
    }
    int_array_free(neighbors);
  }
  polymec_free(R_par);
  polymec_free(x_par);

  // Create a neighbor pairing.
  int num_pairs = (int)pair_array->size/2;
  neighbor_pairing_t* neighbors =
    neighbor_pairing_new("Distance-based point pairs",
                         num_pairs, pair_array->data, ex);

  // Set the number of ghost points referred to within the neighbor pairing.
  *num_ghost_points = (int)(kd_tree_size(tree) - points->num_points);

  // Clean up.
  int_array_release_data_and_free(pair_array); // Release control of data.
  kd_tree_free(tree);

  return neighbors;
}
