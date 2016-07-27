// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "core/block_diagonal_matrix.h"
#include "core/sparse_local_matrix.h"
#include "core/dense_local_matrix.h"
#include "integrators/cpr_differencer.h"

static int F_ident(void* bs_p, real_t t, real_t* x, real_t* F)
{
  int block_size = *((int*)bs_p);
  ASSERT(block_size > 0);
  memcpy(F, x, sizeof(real_t) * block_size * 10);
  return 0;
}

void test_identity_jacobian_with_bs(void** state, int block_size)
{
  adj_graph_t* sparsity = adj_graph_new(MPI_COMM_SELF, 10);
  adj_graph_t* block_sparsity = adj_graph_new_with_block_size(sparsity, block_size);
  // Note: block_sparsity is consumed by diff!
  cpr_differencer_t* diff = cpr_differencer_new(MPI_COMM_SELF,
                                                &block_size,
                                                F_ident,
                                                NULL, 
                                                NULL,
                                                block_sparsity, 
                                                10*block_size, 0);

  // x vector.
  real_t x[block_size * 10];
  for (int i = 0; i < block_size * 10; ++i)
    x[i] = 1.0*i;

  // Try (1) a block diagonal matrix and (2) a sparse matrix.
  local_matrix_t* matrices[2];
  matrices[0] = block_diagonal_matrix_new(10, block_size);
  matrices[1] = sparse_local_matrix_new(block_sparsity);

  for (int m = 0; m < 2; ++m)
  {
    local_matrix_t* J = matrices[m];
    cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, J);
    local_matrix_fprintf(J, stdout); 

    for (int i = 0; i < block_size * 10; ++i)
    {
      for (int j = 0; j < block_size * 10; ++j)
      {
        if (j == i)
        {
          assert_true(local_matrix_value(J, i, j) == 1.0);
        }
        else
        {
          assert_true(local_matrix_value(J, i, j) == 0.0);
        }
      }
    }
    local_matrix_free(J);
  }

  // Clean up.
  adj_graph_free(sparsity);
  cpr_differencer_free(diff);
}

void test_identity_jacobian(void** state)
{
  test_identity_jacobian_with_bs(state, 1);
  test_identity_jacobian_with_bs(state, 2);
  test_identity_jacobian_with_bs(state, 3);
  test_identity_jacobian_with_bs(state, 4);
}

static int F_dense(void* bs_p, real_t t, real_t* x, real_t* F)
{
  int block_size = *((int*)bs_p);
  ASSERT(block_size > 0);
  memset(F, 0, sizeof(real_t) * block_size * 10);
  for (int i = 0; i < block_size * 10; ++i)
    for (int j = 0; j < block_size * 10; ++j)
      F[i] += x[j]; 
  return 0;
}

void test_dense_jacobian_with_bs(void** state, int block_size)
{
  // Make a dense sparsity graph.
  adj_graph_t* sparsity = dense_adj_graph_new(MPI_COMM_SELF, 10, 0);
  adj_graph_t* block_sparsity = adj_graph_new_with_block_size(sparsity, block_size);
  adj_graph_free(sparsity);
  // Note: block_sparsity is consumed by diff!
  cpr_differencer_t* diff = cpr_differencer_new(MPI_COMM_SELF,
                                                &block_size,
                                                F_dense,
                                                NULL,
                                                NULL,
                                                block_sparsity, 
                                                10*block_size, 0);

  // x vector.
  real_t x[block_size * 10];
  for (int i = 0; i < block_size * 10; ++i)
    x[i] = 1.0*i;

  // Try a block diagonal matrix.
  local_matrix_t* D = block_diagonal_matrix_new(10, block_size);
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, D);
  local_matrix_fprintf(D, stdout); 
  for (int i = 0; i < block_size * 10; ++i)
  {
    int block_row = i / block_size;
    for (int j = 0; j < block_size * 10; ++j)
    {
      int block_col = j / block_size;
      if (block_col == block_row)
      {
        assert_true(local_matrix_value(D, i, j) == 1.0);
      }
      else
      {
        assert_true(local_matrix_value(D, i, j) == 0.0);
      }
    }
  }
  local_matrix_free(D);

  // Now try sparse and dense matrices.
  local_matrix_t* A = sparse_local_matrix_new(block_sparsity);
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, A);
  local_matrix_fprintf(A, stdout);
  for (int i = 0; i < block_size * 10; ++i)
  {
    for (int j = 0; j < block_size * 10; ++j)
    {
      assert_true(local_matrix_value(A, i, j) == 1.0);
    }
  }
  local_matrix_free(A);

  A = dense_local_matrix_new(10*block_size);
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, A);
  local_matrix_fprintf(A, stdout);
  for (int i = 0; i < block_size * 10; ++i)
  {
    for (int j = 0; j < block_size * 10; ++j)
    {
      assert_true(local_matrix_value(A, i, j) == 1.0);
    }
  }
  local_matrix_free(A);

  // Clean up.
  cpr_differencer_free(diff);
}

void test_dense_jacobian(void** state)
{
  test_dense_jacobian_with_bs(state, 1);
  test_dense_jacobian_with_bs(state, 2);
  test_dense_jacobian_with_bs(state, 3);
  test_dense_jacobian_with_bs(state, 4);
}

static int F_asym(void* null, real_t t, real_t* x, real_t* F)
{
  memset(F, 0, sizeof(real_t) * 10);
  int counter = 1;
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j, ++counter)
      F[i] += counter * x[j]; 
  return 0;
}

void test_asymmetric_jacobian(void** state)
{
  // Make a dense sparsity graph.
  adj_graph_t* sparsity = dense_adj_graph_new(MPI_COMM_SELF, 10, 0);
  // Note: sparsity is consumed by diff!
  cpr_differencer_t* diff = cpr_differencer_new(MPI_COMM_SELF,
                                                NULL,
                                                F_asym,
                                                NULL,
                                                NULL,
                                                sparsity, 
                                                10, 0);

  // x vector.
  real_t x[10];
  for (int i = 0; i < 10; ++i)
    x[i] = 1.0*i;

  // Try a block diagonal matrix.
  local_matrix_t* D = block_diagonal_matrix_new(10, 1);
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, D);
  local_matrix_fprintf(D, stdout); 
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      if (i == j)
      {
        assert_true(local_matrix_value(D, i, j) == (10.0*i + 1.0*j + 1.0));
      }
      else
      {
        assert_true(local_matrix_value(D, i, j) == 0.0);
      }
    }
  }
  local_matrix_free(D);

  // Now try sparse and dense matrices.
  local_matrix_t* A = sparse_local_matrix_new(sparsity);
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, A);
  local_matrix_fprintf(A, stdout);
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      assert_true(local_matrix_value(A, i, j) == (10.0*i + 1.0*j + 1.0));
    }
  }
  local_matrix_free(A);

  A = dense_local_matrix_new(10);
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, A);
  local_matrix_fprintf(A, stdout);
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {
      assert_true(local_matrix_value(A, i, j) == (10.0*i + 1.0*j + 1.0));
    }
  }
  local_matrix_free(A);

  // Clean up.
  cpr_differencer_free(diff);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_identity_jacobian),
    cmocka_unit_test(test_dense_jacobian),
    cmocka_unit_test(test_asymmetric_jacobian)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
