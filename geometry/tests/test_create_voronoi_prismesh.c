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
#include "geometry/create_voronoi_prismesh.h"

static void test_create_voronoi_prismesh_in_bbox(void** state)
{
  // Generate 500 uniformly distributed random points in the plane.
  size_t N = 500;
  rng_t* rng = host_rng_new();
  point2_t generators[N];
  for (size_t i = 0; i < N; ++i)
    generators[i].x = rng_uniform(rng);

  // Here's a bounding box to contain the mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};

  // Generate a local "NxN" voronoi prismesh from these points within a bounding box, and verify it.
  {
    prismesh_t* mesh = create_voronoi_prismesh_in_bbox(MPI_COMM_SELF, &bbox, generators, N, N);

    assert_int_equal(1, prismesh_num_layers(mesh));
    assert_int_equal(N, prismesh_num_columns(mesh));
    assert_int_equal(N, prismesh_num_vertical_cells(mesh));

    // Leaf through the polygons in the mesh.
    polygon_t* polygons[N];
    for (size_t i = 0; i < N; ++i)
    {
      polygon_t* poly = prismesh_polygon(mesh, i);
      assert_true(poly != NULL);
      polygons[i] = poly;
    }

    // Loop over the layers in the mesh.
    int pos = 0;
    size_t layers_visited = 0;
    prismesh_layer_t* layer;
    while (prismesh_next_layer(mesh, &pos, &layer))
    {
      assert_int_equal(N, layer->columns->size);

      // Check layer lower/upper boundaries.
      assert_true(reals_equal(layer->z1, 0.0));
      assert_true(reals_equal(layer->z2, 1.0));

      // Visit all the columns.
      int ppos = 0;
      size_t cols_visited = 0;
      prismesh_column_t* column;
      while (prismesh_layer_next_column(layer, &ppos, &column))
      {
        assert_true(polygons[column->index] == column->polygon);
        assert_int_equal(N, column->num_cells);
        assert_true(column->neighbors != NULL);
        size_t ne = polygon_num_edges(column->polygon);
        for (size_t e = 0; e < ne; ++e)
          assert_true(column->neighbors[e] != NULL);
        ++cols_visited;
      }
      assert_int_equal(N, cols_visited);

      ++layers_visited;
    }
    assert_int_equal(1, layers_visited);

    // Clean up.
    prismesh_free(mesh);
  }
  
  // Generate a distributed "NxN" voronoi prismesh from these points within a bounding box, and verify it.
  {
    prismesh_t* mesh = create_voronoi_prismesh_in_bbox(MPI_COMM_WORLD, &bbox, generators, N, N);

    size_t num_layers = prismesh_num_layers(mesh);
    assert_true(num_layers >= 1);
    assert_int_equal(N, prismesh_num_columns(mesh));
    assert_int_equal(N, prismesh_num_vertical_cells(mesh));

    // Leaf through the polygons in the mesh.
    polygon_t* polygons[N];
    for (size_t i = 0; i < N; ++i)
    {
      polygon_t* poly = prismesh_polygon(mesh, i);
      assert_true(poly != NULL);
      polygons[i] = poly;
    }

    // Loop over the layers in the mesh.
    int pos = 0;
    size_t layers_visited = 0;
    prismesh_layer_t* layer;
    while (prismesh_next_layer(mesh, &pos, &layer))
    {
      assert_int_equal(N, layer->columns->size);

      // Check layer lower/upper boundaries.
//      assert_true(reals_equal(layer->z1, 0.0));
//      assert_true(reals_equal(layer->z2, 1.0));

      // Visit all the columns.
      int ppos = 0;
      size_t cols_visited = 0;
      prismesh_column_t* column;
      while (prismesh_layer_next_column(layer, &ppos, &column))
      {
        assert_true(polygons[column->index] == column->polygon);
        assert_int_equal(N, column->num_cells);
        assert_true(column->neighbors != NULL);
        size_t ne = polygon_num_edges(column->polygon);
        for (size_t e = 0; e < ne; ++e)
          assert_true(column->neighbors[e] != NULL);
        ++cols_visited;
      }
      assert_int_equal(N, cols_visited);

      ++layers_visited;
    }
    assert_int_equal(num_layers, layers_visited);

    // Clean up.
    prismesh_free(mesh);
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_voronoi_prismesh_in_bbox),
    cmocka_unit_test(test_create_voronoi_prismesh_in_polygon)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
