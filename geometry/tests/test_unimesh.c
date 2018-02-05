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
#include "geometry/unimesh.h"

static void test_ctors(void** state) 
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};

  unimesh_t* mesh = unimesh_new(MPI_COMM_SELF, &bbox, 4, 4, 4, 10, 10, 10,
                                false, false, false);
  assert_true(unimesh_comm(mesh) == MPI_COMM_SELF);
  assert_int_equal(4*4*4, unimesh_num_patches(mesh));

  real_t dx, dy, dz;
  unimesh_get_spacings(mesh, &dx, &dy, &dz);
  assert_true(reals_equal(dx, 1.0/40));
  assert_true(reals_equal(dy, 1.0/40));
  assert_true(reals_equal(dz, 1.0/40));

  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  assert_int_equal(4, npx);
  assert_int_equal(4, npy);
  assert_int_equal(4, npz);

  int nx, ny, nz;
  unimesh_get_patch_size(mesh, &nx, &ny, &nz);
  assert_int_equal(10, nx);
  assert_int_equal(10, ny);
  assert_int_equal(10, nz);

  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  assert_false(x_periodic);
  assert_false(y_periodic);
  assert_false(z_periodic);

  unimesh_free(mesh);

  mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 4, 4, 4, 10, 10, 10, 
                     true, true, true);
  assert_true(unimesh_comm(mesh) == MPI_COMM_WORLD);
  assert_true(unimesh_num_patches(mesh) <= 4*4*4);

  unimesh_get_spacings(mesh, &dx, &dy, &dz);
  assert_true(reals_equal(dx, 1.0/40));
  assert_true(reals_equal(dy, 1.0/40));
  assert_true(reals_equal(dz, 1.0/40));

  unimesh_get_extents(mesh, &npx, &npy, &npz);
  assert_int_equal(4, npx);
  assert_int_equal(4, npy);
  assert_int_equal(4, npz);

  unimesh_get_patch_size(mesh, &nx, &ny, &nz);
  assert_int_equal(10, nx);
  assert_int_equal(10, ny);
  assert_int_equal(10, nz);

  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  assert_true(x_periodic);
  assert_true(y_periodic);
  assert_true(z_periodic);

  unimesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctors)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
