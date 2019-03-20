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
#include "geometry/create_hex_planar_polymesh.h"

extern void write_plot_planar_polymesh_py(planar_polymesh_t* mesh,
                                          const char* script_name);

static void test_create_hex_planar_polymesh(void** state,
                                            int radius,
                                            int num_cells,
                                            int num_edges,
                                            int num_nodes)
{
  real_t h = 0.1;
  planar_polymesh_t* mesh = create_hex_planar_polymesh(radius, h);

  // Check its connectivity.
  assert_true(mesh->num_cells == num_cells);
  assert_true(mesh->num_edges == num_edges);
  assert_true(mesh->num_nodes == num_nodes);

  // Write a Python/Matplotlib script to plot the mesh.
  char script[FILENAME_MAX+1];
  snprintf(script, FILENAME_MAX, "plot_r%d_planar_hexmesh.py", radius);
  write_plot_planar_polymesh_py(mesh, script);

  // Verify its topology.
  assert_true(planar_polymesh_is_valid(mesh, NULL));

  planar_polymesh_free(mesh);
}

static void test_create_single_cell_hex_planar_polymesh(void** state)
{
  test_create_hex_planar_polymesh(state, 0, 1, 6, 6);
}

static void test_create_r1_hex_planar_polymesh(void** state)
{
  test_create_hex_planar_polymesh(state, 1, 7, 30, 27);
}

static void test_create_r5_hex_planar_polymesh(void** state)
{
  test_create_hex_planar_polymesh(state, 5, 91, 306, 232);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_create_single_cell_hex_planar_polymesh),
    cmocka_unit_test(test_create_r1_hex_planar_polymesh),
    cmocka_unit_test(test_create_r5_hex_planar_polymesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
