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
#include "core/graph_from_mesh_cells.h"
#include "core/create_uniform_mesh.h"

void test_graph_from_uniform_mesh_cells(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  bbox_t box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* m = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &box);
  adj_graph_t* g = graph_from_mesh_cells(m);
  int first = adj_graph_first_vertex(g);
  int last = adj_graph_last_vertex(g);
  int interior_cells = 0, plane_cells = 0, edge_cells = 0, corner_cells = 0;
  for (int v = first; v <= last; ++v)
  {
    int ne = adj_graph_num_edges(g, v);
    if (ne == 6) ++interior_cells;
    else if (ne == 5) ++plane_cells;
    else if (ne == 4) ++edge_cells;
    else 
    {
      assert_int_equal(3, ne);
      ++corner_cells;
    }
  }
  assert_int_equal(8, corner_cells);
  assert_int_equal(8*12, edge_cells);
  assert_int_equal(64*6, plane_cells);
  assert_int_equal(8*8*8, interior_cells);
  adj_graph_free(g);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_graph_from_uniform_mesh_cells)
  };
  return run_tests(tests);
}
