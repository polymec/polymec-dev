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

#include "core/tuple.h"
#include "geometry/blockmesh.h"

//#include "geometry/tests/create_cubed_sphere.h"
#include "geometry/tests/create_multiblock_mesh.h"

static void test_serial_ctor(void** state)
{
  blockmesh_t* mesh = create_multiblock_mesh(MPI_COMM_SELF,
                                             2, 2, 10, 10, 0.9, 1.0);
//  blockmesh_t* mesh = create_cubed_sphere(MPI_COMM_SELF,
//                                          2, 2, 10, 10, 0.9, 1.0);
  assert_true(blockmesh_comm(mesh) == MPI_COMM_SELF);
  assert_int_equal(4, blockmesh_num_blocks(mesh));
  for (int b = 0; b < 4; ++b)
  {
    assert_true(blockmesh_block_is_connected(mesh, b, UNIMESH_X1_BOUNDARY));
    assert_true(blockmesh_block_is_connected(mesh, b, UNIMESH_X2_BOUNDARY));
  }
  blockmesh_free(mesh);
}

static void test_parallel_ctor(void** state)
{
  blockmesh_t* mesh = create_multiblock_mesh(MPI_COMM_WORLD,
                                             2, 2, 10, 10, 0.9, 1.0);
  assert_true(blockmesh_comm(mesh) == MPI_COMM_WORLD);
  assert_int_equal(4, blockmesh_num_blocks(mesh));
  for (int b = 0; b < 4; ++b)
  {
    assert_true(blockmesh_block_is_connected(mesh, b, UNIMESH_X1_BOUNDARY));
    assert_true(blockmesh_block_is_connected(mesh, b, UNIMESH_X2_BOUNDARY));
  }
  blockmesh_free(mesh);
}

static void test_next_block(void** state)
{
  blockmesh_t* mesh = create_multiblock_mesh(MPI_COMM_SELF,
                                             2, 2, 10, 10, 0.9, 1.0);
  int pos = 0, block_index;
  unimesh_t* block;
  while (blockmesh_next_block(mesh, &pos, &block_index, &block))
    assert_true(block == blockmesh_block(mesh, block_index));
  blockmesh_free(mesh);
}

static void test_repartition(void** state)
{
  blockmesh_t* mesh = create_multiblock_mesh(MPI_COMM_SELF,
                                             2, 2, 10, 10, 0.9, 1.0);
  repartition_blockmesh(&mesh, NULL, 0.05, NULL, 0);
  blockmesh_free(mesh);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_serial_ctor),
    cmocka_unit_test(test_parallel_ctor),
    cmocka_unit_test(test_next_block),
    cmocka_unit_test(test_repartition)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
