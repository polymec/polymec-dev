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
#include "geometry/polymesh.h"
#include "geometry/create_uniform_polymesh.h"
#include "solvers/matrix_sparsity.h"

static void test_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  index_t row_dist[nprocs+1];
  row_dist[0] = 0;
  for (int i = 0; i < nprocs; ++i)
    row_dist[i+1] = row_dist[i] + 100;
  matrix_sparsity_t* S = matrix_sparsity_new(MPI_COMM_WORLD, row_dist);
  assert_int_equal(100*nprocs, matrix_sparsity_num_global_rows(S));
  assert_true(matrix_sparsity_num_local_rows(S) <= matrix_sparsity_num_global_rows(S));
  index_t* dist = matrix_sparsity_row_distribution(S);
  for (int i = 0; i <= nprocs; ++i)
  {
    assert_int_equal(row_dist[i], dist[i]);
  }
  matrix_sparsity_fprintf(S, stdout);

  matrix_sparsity_t* S1 = matrix_sparsity_clone(S);
  assert_int_equal(100*nprocs, matrix_sparsity_num_global_rows(S1));

  matrix_sparsity_free(S1);
  matrix_sparsity_free(S);
}

static void test_graph_constructor(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, 10, 10, 10, &bbox);
  adj_graph_t* g = graph_from_polymesh_cells(mesh);
  exchanger_t* ex = polymesh_exchanger(mesh);
  matrix_sparsity_t* S = matrix_sparsity_from_graph(g, ex);
  assert_true(matrix_sparsity_comm(S) == mesh->comm);
  assert_int_equal(1000, matrix_sparsity_num_global_rows(S));

  // Make a sparsity with 2x2 blocks.
  matrix_sparsity_t* S1 = matrix_sparsity_with_block_size(S, 2);
  assert_int_equal(2000, matrix_sparsity_num_global_rows(S1));

  // Make a sparsity with variable-sized blocks.
  size_t block_sizes[1000];
  for (int i = 0; i < 1000; ++i)
    block_sizes[i] = i % 4 + 1;
  matrix_sparsity_t* S11 = matrix_sparsity_with_block_sizes(S, block_sizes);

  matrix_sparsity_free(S11);
  matrix_sparsity_free(S1);
  matrix_sparsity_free(S);
  adj_graph_free(g);
  polymesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_constructor),
    cmocka_unit_test(test_graph_constructor)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
