// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/matrix_sparsity.h"

struct matrix_sparsity_t 
{
  MPI_Comm comm;
  int rank, nproc;
  index_t* row_dist; // Global row distribution.
  size_t num_local_rows; // Number of locally stored rows.
  size_t num_global_rows; // Number of locally stored rows.
  index_t* columns; // Column indices in compressed row storage format.
  size_t* offsets; // Column offsets array.
  size_t columns_cap; // Current capacity of columns array.
};

matrix_sparsity_t* matrix_sparsity_new(MPI_Comm comm, 
                                       index_t* row_distribution)
{
  matrix_sparsity_t* sparsity = polymec_malloc(sizeof(matrix_sparsity_t));
  sparsity->comm = comm;
  MPI_Comm_size(comm, &sparsity->nproc);
  MPI_Comm_rank(comm, &sparsity->rank);
  sparsity->row_dist = polymec_malloc(sizeof(index_t) * (sparsity->nproc+1));
  memcpy(sparsity->row_dist, row_distribution, sizeof(index_t) * (sparsity->nproc+1));
  sparsity->num_local_rows = row_distribution[sparsity->rank+1] - row_distribution[sparsity->rank];
  sparsity->num_global_rows = row_distribution[sparsity->nproc];

  // We provide diagonals "for free."
  int init_row_cap = 4; // Initial allocation for columns per row (incl. diagonal)
  sparsity->columns_cap = init_row_cap * sparsity->num_local_rows;
  sparsity->columns = polymec_calloc(sparsity->columns_cap, sizeof(index_t));
  sparsity->offsets = polymec_malloc(sizeof(size_t) * (sparsity->num_local_rows + 1));
  sparsity->offsets[0] = 0;
  for (size_t i = 0; i < sparsity->num_local_rows; ++i)
  {
    sparsity->columns[init_row_cap*i] = i;
    sparsity->offsets[i+1] = sparsity->offsets[i] + init_row_cap;
  }
  return sparsity;
}

matrix_sparsity_t* matrix_sparsity_with_block_size(matrix_sparsity_t* sparsity,
                                                   size_t block_size)
{
  ASSERT(block_size > 0);
  int nr = (int)sparsity->num_local_rows;
  size_t block_sizes[nr];
  for (size_t i = 0; i < nr; ++i)
    block_sizes[i] = block_size;
  return matrix_sparsity_with_block_sizes(sparsity, block_sizes);
}

matrix_sparsity_t* matrix_sparsity_with_block_sizes(matrix_sparsity_t* sparsity,
                                                    size_t* block_sizes)
{
  // We will create a new matrix sparsity pattern with the same layout as 
  // the original one.
  matrix_sparsity_t* block_sp = polymec_malloc(sizeof(matrix_sparsity_t));
  block_sp->comm = sparsity->comm;
  block_sp->nproc = sparsity->nproc;
  block_sp->rank = sparsity->rank;
  block_sp->row_dist = polymec_malloc(sizeof(index_t) * (block_sp->nproc+1));

  // Count up local and global rows, and distribute them accordingly.
  size_t my_num_rows = 0;
  size_t num_rows[block_sp->nproc];
  for (size_t i = 0; i < sparsity->num_local_rows; ++i)
    my_num_rows += block_sizes[i];
  MPI_Allgather(&my_num_rows, 1, MPI_SIZE_T, 
                num_rows, 1, MPI_SIZE_T, 
                block_sp->comm);
  block_sp->row_dist[0] = 0;
  for (int p = 0; p < block_sp->nproc; ++p)
    block_sp->row_dist[p+1] = block_sp->row_dist[p] + num_rows[p];
  block_sp->num_local_rows = block_sp->row_dist[block_sp->rank+1] - block_sp->row_dist[block_sp->rank];
  block_sp->num_global_rows = block_sp->row_dist[block_sp->nproc];

  // Allocate memory for columns.
  block_sp->columns_cap = 0;
  for (size_t i = 0; i < sparsity->num_local_rows; ++i)
  {
    size_t bs = block_sizes[i];
    block_sp->columns_cap += bs * bs * (sparsity->offsets[i+1] - sparsity->offsets[i]);
  }
  block_sp->columns = polymec_malloc(sizeof(index_t) * block_sp->columns_cap);
  block_sp->offsets = polymec_malloc(sizeof(size_t) * (block_sp->num_local_rows + 1));

  // Assign column indices.
  block_sp->offsets[0] = 0;
  index_t row = 0;
  for (size_t i = 0; i < sparsity->num_local_rows; ++i)
  {
    size_t bs = block_sizes[i];
    for (size_t ii = 0; ii < bs; ++ii, ++row)
    {
      size_t row_offset = block_sp->offsets[row];
      size_t l = 0;
      for (index_t j = sparsity->offsets[i]; j < sparsity->offsets[i+1]; ++j)
      {
        index_t s_col = sparsity->columns[j];
        for (int jj = 0; jj < bs; ++jj, ++l)
          block_sp->columns[row_offset+l] = bs * s_col + jj;
      }
      block_sp->offsets[row+1] = block_sp->offsets[row] + l;
    }
  }
  ASSERT(row == block_sp->num_local_rows);
  return block_sp;
}

matrix_sparsity_t* matrix_sparsity_from_graph(adj_graph_t* graph,
                                              exchanger_t* ex)
{
  MPI_Comm comm = adj_graph_comm(graph);
  ASSERT((ex != NULL) || (comm == MPI_COMM_SELF));
  index_t* row_dist = adj_graph_vertex_dist(graph);
  matrix_sparsity_t* sparsity = matrix_sparsity_new(comm, row_dist);
  size_t num_local_rows = adj_graph_num_vertices(graph);

  // Build a mapping from ghost indices to global indices.
  int_index_unordered_map_t* globals = int_index_unordered_map_new();
  if (sparsity->nproc > 1)
  {
    size_t num_padded_rows = (size_t)(adj_graph_max_vertex_index(graph) + 1);
    index_t* indices = polymec_malloc(sizeof(index_t) * num_padded_rows);
    for (int i = 0; i < (int)num_local_rows; ++i)
      indices[i] = row_dist[sparsity->rank] + i;
    exchanger_exchange(ex, indices, 1, 0, MPI_INDEX_T);

    for (int i = (int)num_local_rows; i < (int)num_padded_rows; ++i)
      int_index_unordered_map_insert(globals, i, indices[i]);
    polymec_free(indices);
  }

  // Create off-diagonal columns for each of the graph's edges.
  int rpos = 0;
  index_t row, first_row = sparsity->row_dist[sparsity->rank];
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    // Set the number of columns.
    int v = (int)(row - first_row);
    size_t ne = adj_graph_num_edges(graph, v);
    matrix_sparsity_set_num_columns(sparsity, row, ne + 1);

    // Now add the edges.
    index_t* columns = matrix_sparsity_columns(sparsity, row);
    columns[0] = row;
    int epos = 0, e = 1, v1;
    while (adj_graph_next_edge(graph, v, &epos, &v1))
    {
      if (v1 >= num_local_rows)
        columns[e] = *int_index_unordered_map_get(globals, v1);
      else
        columns[e] = row_dist[sparsity->rank] + v1;
      ++e;
    }
  }

  // Clean up.
  int_index_unordered_map_free(globals);

  return sparsity;
}

matrix_sparsity_t* matrix_sparsity_clone(matrix_sparsity_t* sparsity)
{
  matrix_sparsity_t* clone = polymec_malloc(sizeof(matrix_sparsity_t));
  clone->comm = sparsity->comm;
  clone->nproc = sparsity->nproc;
  clone->rank = sparsity->rank;
  clone->row_dist = polymec_malloc(sizeof(index_t) * (sparsity->nproc+1));
  memcpy(clone->row_dist, sparsity->row_dist, sizeof(index_t) * (sparsity->nproc+1));
  clone->num_local_rows = sparsity->num_local_rows;
  clone->num_global_rows = sparsity->num_global_rows;
  clone->offsets = polymec_malloc(sizeof(size_t) * (clone->num_local_rows + 1));
  memcpy(clone->offsets, sparsity->offsets, sizeof(size_t) * (clone->num_local_rows + 1));
  clone->columns_cap = sparsity->columns_cap;
  clone->columns = polymec_malloc(sizeof(index_t) * clone->columns_cap);
  memcpy(clone->columns, sparsity->columns, sizeof(index_t) * (clone->offsets[clone->num_local_rows]));
  return clone;
}

void matrix_sparsity_free(matrix_sparsity_t* sparsity)
{
  polymec_free(sparsity->offsets);
  polymec_free(sparsity->columns);
  polymec_free(sparsity->row_dist);
  polymec_free(sparsity);
}

MPI_Comm matrix_sparsity_comm(matrix_sparsity_t* sparsity)
{
  return sparsity->comm;
}

size_t matrix_sparsity_num_global_rows(matrix_sparsity_t* sparsity)
{
  return sparsity->num_global_rows;
}

size_t matrix_sparsity_num_local_rows(matrix_sparsity_t* sparsity)
{
  return sparsity->num_local_rows;
}

index_t* matrix_sparsity_row_distribution(matrix_sparsity_t* sparsity)
{
  return sparsity->row_dist;
}

size_t matrix_sparsity_num_nonzeros(matrix_sparsity_t* sparsity)
{
  return sparsity->offsets[sparsity->num_local_rows];
}

void matrix_sparsity_set_num_columns(matrix_sparsity_t* sparsity, 
                                     index_t row, 
                                     size_t num_columns)
{
  ASSERT(row < sparsity->num_global_rows);
  index_t local_row = row - sparsity->row_dist[sparsity->rank];
  size_t old_num_columns = sparsity->offsets[local_row+1] - sparsity->offsets[local_row];
  size_t tot_num_columns = sparsity->offsets[sparsity->num_local_rows];
  if (num_columns < old_num_columns)
  {
    size_t num_columns_removed = old_num_columns - num_columns;
    for (size_t i = sparsity->offsets[local_row]; i < tot_num_columns - num_columns_removed; ++i)
      sparsity->columns[i] = sparsity->columns[i + num_columns_removed];
    for (index_t i = local_row + 1; i <= (index_t)sparsity->num_local_rows; ++i)
      sparsity->offsets[i] -= num_columns_removed;
  }
  else if (num_columns > old_num_columns)
  {
    size_t num_columns_added = num_columns - old_num_columns;
    if (tot_num_columns + num_columns_added > sparsity->columns_cap)
    {
      while (tot_num_columns + num_columns_added > sparsity->columns_cap)
        sparsity->columns_cap *= 2;
      sparsity->columns = polymec_realloc(sparsity->columns, sizeof(index_t) * sparsity->columns_cap);
    }
    for (size_t i = tot_num_columns-1; (i >= sparsity->offsets[local_row]) && (i < tot_num_columns); --i)
      sparsity->columns[i + num_columns_added] = sparsity->columns[i];
    for (index_t i = local_row + 1; i <= (index_t)sparsity->num_local_rows; ++i)
      sparsity->offsets[i] += num_columns_added;
  }
}

size_t matrix_sparsity_num_columns(matrix_sparsity_t* sparsity, 
                                   index_t row)
{
  index_t local_row = row - sparsity->row_dist[sparsity->rank];
  return sparsity->offsets[local_row+1] - sparsity->offsets[local_row];
}

index_t* matrix_sparsity_columns(matrix_sparsity_t* sparsity, index_t row)
{
  index_t local_row = row - sparsity->row_dist[sparsity->rank];
  return &sparsity->columns[sparsity->offsets[local_row]];
}

bool matrix_sparsity_contains(matrix_sparsity_t* sparsity, 
                              index_t row, 
                              index_t column)
{
  ASSERT(row < matrix_sparsity_num_global_rows(sparsity));
  ASSERT(column < matrix_sparsity_num_global_rows(sparsity));

  index_t local_row = row - sparsity->row_dist[sparsity->rank];
  for (index_t c = sparsity->offsets[local_row]; c < sparsity->offsets[local_row+1]; ++c)
  {
    if (sparsity->columns[c] == column)
      return true;
  }
  return false;
}

bool matrix_sparsity_next_row(matrix_sparsity_t* sparsity, 
                              int* pos, 
                              index_t* row)
{
  if (*pos >= sparsity->num_local_rows)
    return false;
  *row = sparsity->row_dist[sparsity->rank] + *pos;
  ++(*pos);
  return true;
}

bool matrix_sparsity_next_column(matrix_sparsity_t* sparsity, 
                                 index_t row,
                                 int* pos, 
                                 index_t* column)
{
  index_t local_row = row - sparsity->row_dist[sparsity->rank];
  if ((sparsity->offsets[local_row] + *pos) >= sparsity->offsets[local_row+1])
    return false;
  *column = sparsity->columns[sparsity->offsets[local_row] + *pos];
  ++(*pos);
  return true;
}

void matrix_sparsity_fprintf(matrix_sparsity_t* sparsity, FILE* stream)
{
  if (stream == NULL) return;
  size_t num_global_rows = matrix_sparsity_num_global_rows(sparsity);
  fprintf(stream, "Matrix sparsity (%d / %d rows locally):\n", 
          (int)sparsity->num_local_rows, (int)num_global_rows);
  index_t begin = sparsity->row_dist[sparsity->rank];
  index_t end = sparsity->row_dist[sparsity->rank+1];
  for (index_t i = begin; i < end; ++i)
  {
    fprintf(stream, " %" PRIu64 ": ", i);
    size_t num_cols = matrix_sparsity_num_columns(sparsity, i);
    index_t* cols = matrix_sparsity_columns(sparsity, i);
    for (size_t j = 0; j < num_cols; ++j)
      fprintf(stream, "%" PRIu64 " ", cols[j]);
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}

// This helper is not part of the "official" Polymec API, but is available for use 
// by library functions. It replaces a global sparsity pattern with a distributed 
// sparsity pattern using the new communicator and row distribution.
void distribute_matrix_sparsity(matrix_sparsity_t** sparsity,
                                MPI_Comm comm,
                                index_t* row_distribution);
void distribute_matrix_sparsity(matrix_sparsity_t** sparsity,
                                MPI_Comm comm,
                                index_t* row_distribution)
{
  matrix_sparsity_t* dist_sparsity = matrix_sparsity_new(comm, row_distribution);

  // Just copy the data from the global pattern to the distributed one.
  int rpos = 0;
  index_t row;
  while (matrix_sparsity_next_row(dist_sparsity, &rpos, &row))
  {
    size_t ncol = matrix_sparsity_num_columns(*sparsity, row);
    matrix_sparsity_set_num_columns(dist_sparsity, row, ncol);
    index_t* columns = matrix_sparsity_columns(dist_sparsity, row);
    int cpos = 0;
    index_t column, c = 0;
    while (matrix_sparsity_next_column(*sparsity, row, &cpos, &column))
    {
      columns[c] = column;
      ++c;
    }
  }

  // Replace sparsity with the distributed one.
  matrix_sparsity_free(*sparsity);
  *sparsity = dist_sparsity;
}

