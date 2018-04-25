// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/polymec.h"
#include "core/adj_graph.h"

static void test_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, 100);
  assert_int_equal(100*nprocs, adj_graph_num_vertices(g));
  adj_graph_free(g);
}

static void test_block_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, 100);
  size_t block_sizes[100*nprocs];
  for (int i = 0; i < 100*nprocs; ++i)
    block_sizes[i] = (size_t)(i % 9 + 1);
  adj_graph_t* bg = adj_graph_new_with_block_sizes(g, block_sizes);
  adj_graph_free(g);
  adj_graph_free(bg);
}

static void test_dense_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  adj_graph_t* g = dense_adj_graph_new(MPI_COMM_WORLD, 10, 10*(nprocs-1));
  assert_int_equal(10*nprocs, adj_graph_num_vertices(g));
  adj_graph_free(g);
}

static void test_distributed_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  index_t vertex_dist[nprocs+1];
  vertex_dist[0] = 0;
  for (int p = 0; p < nprocs; ++p)
    vertex_dist[p+1] = vertex_dist[p] + 1000/nprocs;
  adj_graph_t* g = adj_graph_new_with_dist(MPI_COMM_WORLD, 1000, vertex_dist);
  assert_int_equal(1000/nprocs, adj_graph_num_vertices(g));
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert_int_equal(rank*1000/nprocs, adj_graph_first_vertex(g));
  assert_int_equal((rank+1)*1000/nprocs-1, adj_graph_last_vertex(g));
  adj_graph_free(g);
}

static void test_array_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, 100);
  adj_graph_t* g1 = adj_graph_from_arrays(MPI_COMM_WORLD, 
                                          adj_graph_vertex_dist(g),
                                          adj_graph_adjacency(g),
                                          adj_graph_edge_offsets(g),
                                          false);
  assert_int_equal(100*nprocs, adj_graph_num_vertices(g1));
  adj_graph_free(g);
  adj_graph_free(g1);
}

static void test_sort(void** state)
{
  // Test the sort of a cyclic graph.
  int N = 10;
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, N);
  for (int v = 0; v < N; ++v)
  {
    adj_graph_set_num_edges(g, v, 1);
    int* edges = adj_graph_edges(g, v);
    edges[0] = (v + 1) % N;
  }
  int sorted[N];
  bool result = adj_graph_sort(g, sorted);
  assert_false(result);
  adj_graph_free(g);

  // Now test the sort of an acyclic graph.
  g = adj_graph_new(MPI_COMM_SELF, N);
  for (int v = 0; v < N-1; ++v)
  {
    adj_graph_set_num_edges(g, v, 1);
    int* edges = adj_graph_edges(g, v);
    edges[0] = v + 1;
  }
  result = adj_graph_sort(g, sorted);
  assert_true(result);
  for (int v = 0; v < N; ++v)
  {
    assert_int_equal(sorted[v], v);
  }
  adj_graph_free(g);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_constructor),
    cmocka_unit_test(test_block_constructor),
    cmocka_unit_test(test_distributed_constructor),
    cmocka_unit_test(test_dense_constructor),
    cmocka_unit_test(test_array_constructor),
    cmocka_unit_test(test_sort)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
