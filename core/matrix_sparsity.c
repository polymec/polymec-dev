// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/matrix_sparsity.h"

struct matrix_sparsity_t 
{
  MPI_Comm comm;
  int rank, nproc;
  index_t* row_dist; // Global row distribution.
  index_t num_local_rows; // Number of locally stored rows.
  index_t num_global_rows; // Number of locally stored rows.
  index_t* columns; // Column indices in compressed row storage format.
  index_t* offsets; // Column offsets array.
  index_t columns_cap; // Current capacity of columns array.
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
  sparsity->columns = polymec_malloc(sizeof(index_t) * sparsity->columns_cap);
  sparsity->offsets = polymec_malloc(sizeof(index_t) * (sparsity->num_local_rows + 1));
  sparsity->offsets[0] = 0;
  for (index_t i = 0; i < sparsity->num_local_rows; ++i)
  {
    sparsity->columns[init_row_cap*i] = i;
    sparsity->offsets[i+1] = sparsity->offsets[i] + init_row_cap;
  }
  return sparsity;
}

matrix_sparsity_t* matrix_sparsity_from_graph(adj_graph_t* graph,
                                              exchanger_t* ex)
{
  MPI_Comm comm = adj_graph_comm(graph);
  index_t* row_dist = adj_graph_vertex_dist(graph);
  matrix_sparsity_t* sparsity = matrix_sparsity_new(comm, row_dist);

  // Create off-diagonal columns for each of the graph's edges.
  int rpos = 0, v = 0;
  index_t row;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    // Set the number of columns.
    int ne = adj_graph_num_edges(graph, v);
    matrix_sparsity_set_num_columns(sparsity, row, (index_t)(ne + 1));

    // Now add the edges.
    index_t* columns = matrix_sparsity_columns(sparsity, row);
    int epos = 0, e = 0, v1;
    while (adj_graph_next_edge(graph, v, &epos, &v1))
    {
      columns[e] = row_dist[sparsity->rank] + v1;
      ++e;
    }
  }

  // Now correct the off-process column indices with the exchanger.
  if (sparsity->nproc > 1)
  {
    POLYMEC_NOT_IMPLEMENTED
  }

  return sparsity;
}

matrix_sparsity_t* redistributed_matrix_sparsity(matrix_sparsity_t* sparsity,
                                                 MPI_Comm comm,
                                                 index_t* row_distribution)
{
  matrix_sparsity_t* sparsity1 = matrix_sparsity_new(comm, row_distribution);
  POLYMEC_NOT_IMPLEMENTED
  return sparsity1;
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

index_t matrix_sparsity_num_global_rows(matrix_sparsity_t* sparsity)
{
  return sparsity->num_global_rows;
}

index_t matrix_sparsity_num_local_rows(matrix_sparsity_t* sparsity)
{
  return sparsity->num_local_rows;
}

index_t* matrix_sparsity_row_distribution(matrix_sparsity_t* sparsity)
{
  return sparsity->row_dist;
}

index_t matrix_sparsity_num_nonzeros(matrix_sparsity_t* sparsity)
{
  return sparsity->offsets[sparsity->num_local_rows];
}

void matrix_sparsity_set_num_columns(matrix_sparsity_t* sparsity, 
                                     index_t row, 
                                     index_t num_columns)
{
  ASSERT(row < sparsity->num_global_rows);
  index_t local_row = row - sparsity->row_dist[sparsity->rank];
  index_t old_num_columns = sparsity->offsets[local_row+1] - sparsity->offsets[local_row];
  index_t tot_num_columns = sparsity->offsets[sparsity->num_local_rows];
  if (num_columns < old_num_columns)
  {
    index_t num_columns_removed = old_num_columns - num_columns;
    for (index_t i = sparsity->offsets[local_row]; i < tot_num_columns - num_columns_removed; ++i)
      sparsity->columns[i] = sparsity->columns[i + num_columns_removed];
    for (index_t i = local_row + 1; i <= sparsity->num_local_rows; ++i)
      sparsity->offsets[i] -= num_columns_removed;
  }
  else if (num_columns > old_num_columns)
  {
    index_t num_columns_added = num_columns - old_num_columns;
    if (tot_num_columns + num_columns_added > sparsity->columns_cap)
    {
      while (tot_num_columns + num_columns_added > sparsity->columns_cap)
        sparsity->columns_cap *= 2;
      sparsity->columns = polymec_realloc(sparsity->columns, sizeof(index_t) * sparsity->columns_cap);
    }
    for (index_t i = tot_num_columns-1; i >= sparsity->offsets[local_row]; --i)
      sparsity->columns[i + num_columns_added] = sparsity->columns[i];
    for (index_t i = local_row + 1; i <= sparsity->num_local_rows; ++i)
      sparsity->offsets[i] += num_columns_added;
  }
}

index_t matrix_sparsity_num_columns(matrix_sparsity_t* sparsity, 
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
  index_t num_global_rows = matrix_sparsity_num_global_rows(sparsity);
  fprintf(stream, "Matrix sparsity (%" PRIu64" /%" PRIu64 " rows locally):\n", 
          sparsity->num_local_rows, num_global_rows);
  index_t begin = sparsity->row_dist[sparsity->rank];
  index_t end = sparsity->row_dist[sparsity->rank+1];
  for (index_t i = begin; i < end; ++i)
  {
    fprintf(stream, " %" PRIu64 ": ", i);
    index_t num_cols = matrix_sparsity_num_columns(sparsity, i);
    index_t* cols = matrix_sparsity_columns(sparsity, i);
    for (index_t j = 0; j < num_cols; ++j)
      fprintf(stream, "%" PRIu64 " ", cols[j]);
    fprintf(stream, "\n");
  }
  fprintf(stream, "\n");
}

