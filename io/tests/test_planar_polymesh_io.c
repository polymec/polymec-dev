// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "core/array_utils.h"
#include "geometry/create_quad_planar_polymesh.h"
#include "geometry/create_hex_planar_polymesh.h"
#include "io/silo_file.h"

// Verifies a mesh read from a file against a reference mesh.
static void verify_planar_polymesh(void** state,
                                   planar_polymesh_t* read_mesh,
                                   planar_polymesh_t* ref_mesh)
{
  assert_int_equal(read_mesh->num_cells, ref_mesh->num_cells);
  assert_int_equal(read_mesh->num_edges, ref_mesh->num_edges);
  assert_int_equal(read_mesh->num_nodes, ref_mesh->num_nodes);
  for (size_t c = 0; c <= read_mesh->num_cells; ++c)
    assert_int_equal(read_mesh->cell_edge_offsets[c], ref_mesh->cell_edge_offsets[c]);
  for (size_t e = 0; e < read_mesh->cell_edge_offsets[read_mesh->num_cells]; ++e)
    assert_int_equal(read_mesh->cell_edges[e], ref_mesh->cell_edges[e]);
  for (size_t e = 0; e < read_mesh->num_edges; ++e)
  {
    assert_int_equal(read_mesh->edge_cells[2*e], ref_mesh->edge_cells[2*e]);
    assert_int_equal(read_mesh->edge_cells[2*e+1], ref_mesh->edge_cells[2*e+1]);
    assert_int_equal(read_mesh->edge_nodes[2*e], ref_mesh->edge_nodes[2*e]);
    assert_int_equal(read_mesh->edge_nodes[2*e+1], ref_mesh->edge_nodes[2*e+1]);
  }
  for (size_t n = 0; n < read_mesh->num_nodes; ++n)
  {
    assert_true(reals_equal(0.0, point2_distance(&read_mesh->nodes[n],
                                                 &ref_mesh->nodes[n])));
  }
}

static void test_plot_quad_mesh(void** state)
{
  // Create a 10x10 planar quad mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  planar_polymesh_t* mesh = create_quad_planar_polymesh(10, 10, &bbox, false, false);

  // Plot it.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "quad_10x10", "", 1, 0, 0.0);
  silo_file_write_planar_polymesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Now read the mesh from the file.
  real_t t;
  silo = silo_file_open(MPI_COMM_WORLD, "quad_10x10", "", 0, &t);
  assert_true(reals_equal(t, 0.0));
  assert_true(silo_file_contains_planar_polymesh(silo, "mesh"));
  planar_polymesh_t* mesh1 = silo_file_read_planar_polymesh(silo, "mesh");

  // Check its numbers.
  verify_planar_polymesh(state, mesh1, mesh);

  silo_file_close(silo);

  // Clean up.
  planar_polymesh_free(mesh);
  planar_polymesh_free(mesh1);
}

static void test_plot_hex_mesh(void** state)
{
  // Create a set of uniform hexes.
  int radius = 5;
  real_t h = 0.1;
  planar_polymesh_t* mesh = create_hex_planar_polymesh(radius, h);

  // Plot it.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "hex_r=5", "", 1, 0, 0.0);
  silo_file_write_planar_polymesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.

  // Now read the mesh from the file.
  real_t t;
  silo = silo_file_open(MPI_COMM_WORLD, "hex_r=5", "", 0, &t);
  assert_true(reals_equal(t, 0.0));
  assert_true(silo_file_contains_planar_polymesh(silo, "mesh"));
  planar_polymesh_t* mesh1 = silo_file_read_planar_polymesh(silo, "mesh");

  // Check its numbers.
  verify_planar_polymesh(state, mesh1, mesh);

  silo_file_close(silo);

  // Clean up.
  planar_polymesh_free(mesh);
  planar_polymesh_free(mesh1);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  silo_enable_compression(1); // create compressed files.
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_plot_quad_mesh),
    cmocka_unit_test(test_plot_hex_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
