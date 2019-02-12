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

extern blockmesh_t* create_cubed_sphere(MPI_Comm comm,
                                        int block_nxy, int block_nz,
                                        int patch_nxy, int patch_nz,
                                        real_t R1, real_t R2);

static void test_serial_ctor(void** state) 
{
  blockmesh_t* mesh = create_cubed_sphere(MPI_COMM_SELF,
                                          2, 2, 10, 10, 0.9, 1.0);
  assert_true(blockmesh_comm(mesh) == MPI_COMM_SELF);
  assert_int_equal(6, blockmesh_num_blocks(mesh));
  blockmesh_free(mesh);
}

static void test_parallel_ctor(void** state) 
{
  blockmesh_t* mesh = create_cubed_sphere(MPI_COMM_WORLD,
                                          2, 2, 10, 10, 0.9, 1.0);
  assert_true(blockmesh_comm(mesh) == MPI_COMM_WORLD);
  assert_int_equal(6, blockmesh_num_blocks(mesh));
  blockmesh_free(mesh);
}

static void test_next_block(void** state) 
{
  blockmesh_t* mesh = create_cubed_sphere(MPI_COMM_WORLD,
                                          2, 2, 10, 10, 0.9, 1.0);
  int pos = 0, b = 0;
  unimesh_t* block;
  bbox_t domain;
  coord_mapping_t* coords;
  while (blockmesh_next_block(mesh, &pos, &block, &domain, &coords))
  {
    assert_true(block == blockmesh_block(mesh, b));
    ++b;
  }
  blockmesh_free(mesh);
}

static void test_repartition(void** state) 
{
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
