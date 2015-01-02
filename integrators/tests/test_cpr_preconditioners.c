// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
  preconditioner_t* precond = block_jacobi_preconditioner_from_function("test", NULL, sys_func, NULL, g, N, bs);
  preconditioner_free(precond);
  adj_graph_t* bg = adj_graph_new_with_block_size(bs, g);
  precond = block_jacobi_preconditioner_from_function("test", NULL, sys_func, NULL, bg, N, bs);
  preconditioner_free(precond);
  adj_graph_free(bg);
  adj_graph_free(g);
}

void test_lu_ctor(void** state)
{
  int N = 10;
  int bs = 2;
  adj_graph_t* g = linear_graph(N);
  preconditioner_t* precond = lu_preconditioner_from_function("test", NULL, sys_func, NULL, g, N, bs);
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
