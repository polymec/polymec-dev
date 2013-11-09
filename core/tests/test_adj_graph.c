// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
