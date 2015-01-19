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
#include "slu_ddefs.h"
#include "core/polymec.h"
#include "core/array_utils.h"
#include "integrators/curtis_powell_reed_preconditioners.h"

static adj_graph_t* linear_graph(int N)
{
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, N);
  for (int i = 0; i < N; ++i)
  {
    if ((i > 0) && (i < N-1))
    {
      adj_graph_set_num_edges(g, i, 2);
      int* edges = adj_graph_edges(g, i);
      edges[0] = i-1;
      edges[1] = i+1;
    }
    else if (i == 0)
    {
      adj_graph_set_num_edges(g, i, 1);
      int* edges = adj_graph_edges(g, i);
      edges[0] = 1;
    }
    else
    {
      adj_graph_set_num_edges(g, i, N-1);
      int* edges = adj_graph_edges(g, i);
      edges[0] = N-2;
    }
  }
  return g;
}

static int sys_func(void* context, real_t t, real_t* x, real_t* F)
{
  return 0;
}

void test_block_jacobi_ctor(void** state)
{
  int N = 10;
  int bs = 2;
  adj_graph_t* g = linear_graph(N);
  preconditioner_t* precond = block_jacobi_preconditioner_from_function("test", NULL, sys_func, NULL, g, N, 0, bs);
  preconditioner_free(precond);
  adj_graph_t* bg = adj_graph_new_with_block_size(bs, g);
  precond = block_jacobi_preconditioner_from_function("test", NULL, sys_func, NULL, bg, N, 0, bs);
  preconditioner_free(precond);
  adj_graph_free(bg);
  adj_graph_free(g);
}

void test_lu_ctor(void** state)
{
  int N = 10;
  int bs = 2;
  adj_graph_t* g = linear_graph(N);
  preconditioner_t* precond = lu_preconditioner_from_function("test", NULL, sys_func, NULL, g, N, 0, bs);
  preconditioner_free(precond);
  adj_graph_free(g);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block_jacobi_ctor),
    unit_test(test_lu_ctor)
  };
  return run_tests(tests);
}
