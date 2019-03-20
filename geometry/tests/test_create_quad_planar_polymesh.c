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
#include "geometry/create_quad_planar_polymesh.h"

extern void write_plot_planar_polymesh_py(planar_polymesh_t* mesh,
                                          const char* script_name);

static void test_create_quad_planar_polymesh(void** state)
{
  // Create a 10x10 uniform planar mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  planar_polymesh_t* mesh = create_quad_planar_polymesh(10, 10, &bbox, false, false);

  // Check its numbers.
  assert_int_equal(10*10, mesh->num_cells);
  assert_int_equal(2*10*11, mesh->num_edges);
  assert_int_equal(11*11, mesh->num_nodes);

  // Verify its topology.
  assert_true(planar_polymesh_is_valid(mesh, NULL));

  // Write a Python/Matplotlib script to plot the mesh.
  write_plot_planar_polymesh_py(mesh, "plot_planar_quadmesh.py");

  planar_polymesh_free(mesh);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_create_quad_planar_polymesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
