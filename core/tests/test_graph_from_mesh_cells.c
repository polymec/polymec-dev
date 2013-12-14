// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "core/graph_from_mesh_cells.h"
#include "core/create_uniform_mesh.h"

static adj_graph_t* graph_from_uniform_mesh()
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  bbox_t box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* m = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &box);
  adj_graph_t* g = graph_from_mesh_cells(m);
  mesh_free(m);
  return g;
}

void test_graph_from_uniform_mesh_cells(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();

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

void test_block_graphs_from_uniform_mesh_cells(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();
  int first = adj_graph_first_vertex(g);
  int last = adj_graph_last_vertex(g);

  // Test with block sizes of 1-5.
  for (int b = 1; b <= 5; ++b)
  {
    // Make a block graph from the single one.
    adj_graph_t* bg = adj_graph_new_with_block_size(b, g);
    int nv = adj_graph_num_vertices(bg);
    assert_int_equal(b * adj_graph_num_vertices(g), nv);

    for (int v = first; v <= last; ++v)
    {
      assert_int_equal(b * adj_graph_num_edges(g, v) + (b-1), adj_graph_num_edges(bg, b*v));
      assert_false(adj_graph_contains_edge(bg, b*v, b*v));
      for (int bb = 1; bb < b; ++bb)
      {
        assert_int_equal(b * adj_graph_num_edges(g, v) + (b-1), adj_graph_num_edges(bg, b*v+bb));
        assert_true(adj_graph_contains_edge(bg, b*v, b*v+bb));
      }
      int pos = 0, other_v;
      while (adj_graph_next_edge(g, v, &pos, &other_v))
      {
        for (int bb = 0; bb < b; ++bb)
        {
          assert_true(adj_graph_contains_edge(bg, b*v, b*other_v+bb));
        }
      }
    }

    adj_graph_free(bg);
  }

  adj_graph_free(g);
}

void test_coloring_against_graph(adj_graph_coloring_t* coloring, adj_graph_t* graph)
{
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int color = 0; color < num_colors; ++color)
  {
    int pos = 0, v;
    while (adj_graph_coloring_next_vertex(coloring, color, &pos, &v))
    {
      // This color consists of those vertices not connected to v.
      int pos = 0, v1;
      while (adj_graph_next_edge(graph, v, &pos, &v1))
        assert_false(adj_graph_coloring_has_vertex(coloring, color, v1));
    }
  }
}

void test_smallest_last_graph_coloring_on_uniform_mesh(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();
  adj_graph_coloring_t* c = adj_graph_coloring_new(g, SMALLEST_LAST);
  test_coloring_against_graph(c, g);
  adj_graph_coloring_free(c);
  adj_graph_free(g);
}

void test_largest_first_graph_coloring_on_uniform_mesh(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();
  adj_graph_coloring_t* c = adj_graph_coloring_new(g, LARGEST_FIRST);
  test_coloring_against_graph(c, g);
  adj_graph_coloring_free(c);
  adj_graph_free(g);
}

void test_incidence_degree_graph_coloring_on_uniform_mesh(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();
  adj_graph_coloring_t* c = adj_graph_coloring_new(g, INCIDENCE_DEGREE);
  test_coloring_against_graph(c, g);
  adj_graph_coloring_free(c);
  adj_graph_free(g);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_graph_from_uniform_mesh_cells),
    unit_test(test_block_graphs_from_uniform_mesh_cells),
    unit_test(test_smallest_last_graph_coloring_on_uniform_mesh),
    unit_test(test_largest_first_graph_coloring_on_uniform_mesh)
//    unit_test(test_incidence_degree_graph_coloring_on_uniform_mesh)
  };
  return run_tests(tests);
}
