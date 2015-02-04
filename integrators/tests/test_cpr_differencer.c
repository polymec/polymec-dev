// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmockery.h"
#include "core/polymec.h"
#include "core/block_diagonal_matrix.h"
#include "core/sparse_local_matrix.h"
#include "integrators/cpr_differencer.h"

static int F_ident(void* bs_p, real_t t, real_t* x, real_t* F)
{
  int block_size = *((int*)bs_p);
  ASSERT(block_size > 0);
  memcpy(F, x, sizeof(real_t) * block_size * 100);
  return 0;
}

void test_identity_jacobian_with_bs(void** state, int block_size)
{
  adj_graph_t* sparsity = adj_graph_new(MPI_COMM_SELF, 100);
  cpr_differencer_t* diff = cpr_differencer_new(MPI_COMM_SELF,
                                                F_ident,
                                                NULL,
                                                &block_size,
                                                sparsity, 
                                                100, 0, block_size);
  adj_graph_free(sparsity);
  local_matrix_t* D = block_diagonal_matrix_new(100, block_size);
  real_t x[block_size * 100];
  for (int i = 0; i < block_size * 100; ++i)
    x[i] = 1.0*i;
  cpr_differencer_compute(diff, 0.0, 1.0, 0.0, 0.0, x, NULL, D);
  cpr_differencer_free(diff);
  local_matrix_fprintf(D, stdout);
  local_matrix_free(D);
}

void test_identity_jacobian(void** state)
{
  test_identity_jacobian_with_bs(state, 1);
  test_identity_jacobian_with_bs(state, 2);
  test_identity_jacobian_with_bs(state, 3);
  test_identity_jacobian_with_bs(state, 4);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_identity_jacobian)
  };
  return run_tests(tests);
}
