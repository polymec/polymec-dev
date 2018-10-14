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
#include "geometry/create_hex_planar_polymesh.h"

extern void write_plot_planar_polymesh_py(planar_polymesh_t* mesh,
                                          const char* script_name);

static void test_create_single_cell_hex_planar_polymesh(void** state)
{
  real_t h = 0.1;
  planar_polymesh_t* mesh = create_hex_planar_polymesh(0, h);

  // Check its connectivity.
  assert_true(mesh->num_cells == 1);
  assert_true(mesh->num_edges == 6);
  assert_true(mesh->num_nodes == 6);

  // Verify its topology.
  assert_true(planar_polymesh_verify_topology(mesh, polymec_error));

  // Write a Python/Matplotlib script to plot the mesh.
  write_plot_planar_polymesh_py(mesh, "plot_r0_planar_hexmesh.py");

  planar_polymesh_free(mesh);
}

static void test_create_r1_hex_planar_polymesh(void** state)
{
  // Create a set of 7 uniform hexes (a "flower").
  real_t h = 0.1;
  planar_polymesh_t* mesh = create_hex_planar_polymesh(1, h);

  // Check its connectivity.
  assert_true(mesh->num_cells == 7);
  assert_true(mesh->num_edges == 30);
  assert_true(mesh->num_nodes == 24);

  // Verify its topology.
  assert_true(planar_polymesh_verify_topology(mesh, polymec_error));

  // Write a Python/Matplotlib script to plot the mesh.
  write_plot_planar_polymesh_py(mesh, "plot_r1_planar_hexmesh.py");

  planar_polymesh_free(mesh);
}

static void test_create_r5_hex_planar_polymesh(void** state)
{
  // Create a set of uniform hexes (radius 5).
  real_t h = 0.1;
  planar_polymesh_t* mesh = create_hex_planar_polymesh(5, h);

  // Check its connectivity.
  assert_true(mesh->num_cells == 91);
  assert_true(mesh->num_edges == 306);
  assert_true(mesh->num_nodes == 111);

  // Verify its topology.
  assert_true(planar_polymesh_verify_topology(mesh, polymec_error));

  // Write a Python/Matplotlib script to plot the mesh.
  write_plot_planar_polymesh_py(mesh, "plot_r5_planar_hexmesh.py");

  planar_polymesh_free(mesh);
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
