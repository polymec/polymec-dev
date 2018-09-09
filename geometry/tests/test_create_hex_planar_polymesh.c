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

static void test_create_hex_planar_polymesh(void** state,
                                            hex_lattice_align_t alignment)
{
  // Create a set of uniform hexes with the given alignment.
  size_t radius = 5;
  real_t h = 0.1;
  planar_polymesh_t* mesh = create_hex_planar_polymesh(alignment, radius, h);

  // Check its connectivity.
  assert_true(mesh->num_cells == 10*10);
  assert_true(mesh->num_edges == 4*10*11);
  assert_true(mesh->num_nodes == 11*11);

  planar_polymesh_free(mesh);
}

static void test_create_x_hex_planar_polymesh(void** state)
{
  test_create_hex_planar_polymesh(state, HEX_LATTICE_X_ALIGNED);
}

static void test_create_y_hex_planar_polymesh(void** state)
{
  test_create_hex_planar_polymesh(state, HEX_LATTICE_Y_ALIGNED);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_x_hex_planar_polymesh),
    cmocka_unit_test(test_create_y_hex_planar_polymesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
