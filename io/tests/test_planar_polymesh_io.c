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
#include "core/array_utils.h"
#include "geometry/create_quad_planar_polymesh.h"
#include "io/silo_file.h"

static void test_plot_quad_mesh(void** state)
{
  // Create a 10x10 planar quad mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  planar_polymesh_t* mesh = create_quad_planar_polymesh(10, 10, &bbox, false, false);

  // Plot it.
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "quad_10x10", "", 1, 0, 0.0);
  silo_file_write_planar_polymesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  planar_polymesh_free(mesh);

  // Now read the mesh from the file.
  real_t t;
  silo = silo_file_open(MPI_COMM_WORLD, "rectilinear_4x4x4", "", 0, &t);
  assert_true(reals_equal(t, 0.0));
  assert_true(silo_file_contains_planar_polymesh(silo, "mesh"));
  mesh = silo_file_read_planar_polymesh(silo, "mesh");

  // Check its numbers.
  assert_int_equal(10*10, mesh->num_cells);
  assert_int_equal(2*10*11, mesh->num_edges);
  assert_int_equal(11*11, mesh->num_nodes);

  silo_file_close(silo);

  // Clean up.
  planar_polymesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  silo_enable_compression(1); // create compressed files.
  set_log_level(LOG_DEBUG);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_plot_quad_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
