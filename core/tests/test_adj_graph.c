// Copyright (c) 2012-2014, Jeffrey N. Johnson
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
#include "core/polymec.h"
#include "core/adj_graph.h"

void test_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, 100);
  assert_int_equal(100*nprocs, adj_graph_num_vertices(g));
  adj_graph_free(g);
}

void test_distributed_constructor(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int vertex_dist[nprocs+1];
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

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_constructor),
    unit_test(test_distributed_constructor)
  };
  return run_tests(tests);
}
