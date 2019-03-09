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
#include "geometry/unimesh.h"

// Patch dimensions.
static const int nx = 4;
static const int ny = 4;
static const int nz = 4;

static void test_ctors(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};

  unimesh_t* mesh = unimesh_new(MPI_COMM_SELF, &bbox, 4, 4, 4, nx, ny, nz,
                                false, false, false);
  assert_true(unimesh_comm(mesh) == MPI_COMM_SELF);
  assert_int_equal(4*4*4, unimesh_num_patches(mesh));

  real_t dx, dy, dz;
  unimesh_get_spacings(mesh, &dx, &dy, &dz);
  assert_true(reals_equal(dx, 1.0/(4*nx)));
  assert_true(reals_equal(dy, 1.0/(4*ny)));
  assert_true(reals_equal(dz, 1.0/(4*nz)));

  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  assert_int_equal(4, npx);
  assert_int_equal(4, npy);
  assert_int_equal(4, npz);

  int nx_, ny_, nz_;
  unimesh_get_patch_size(mesh, &nx_, &ny_, &nz_);
  assert_int_equal(nx, nx_);
  assert_int_equal(ny, ny_);
  assert_int_equal(nz, nz_);

  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  assert_false(x_periodic);
  assert_false(y_periodic);
  assert_false(z_periodic);

  unimesh_free(mesh);

  mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 4, 4, 4, nx, ny, nz,
                     true, true, true);
  assert_true(unimesh_comm(mesh) == MPI_COMM_WORLD);
  assert_true(unimesh_num_patches(mesh) <= 4*4*4);

  unimesh_get_spacings(mesh, &dx, &dy, &dz);
  assert_true(reals_equal(dx, 1.0/(4*nx)));
  assert_true(reals_equal(dy, 1.0/(4*ny)));
  assert_true(reals_equal(dz, 1.0/(4*nz)));

  unimesh_get_extents(mesh, &npx, &npy, &npz);
  assert_int_equal(4, npx);
  assert_int_equal(4, npy);
  assert_int_equal(4, npz);

  unimesh_get_patch_size(mesh, &nx_, &ny_, &nz_);
  assert_int_equal(nx, nx_);
  assert_int_equal(ny, ny_);
  assert_int_equal(nz, nz_);

  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  assert_true(x_periodic);
  assert_true(y_periodic);
  assert_true(z_periodic);

  unimesh_free(mesh);
}

// Guarantees that unimesh_next_patch traverses patches in lexicographic order,
// especially in parallel environments.
static void test_next_patch(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};

  unimesh_t* mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 4, 4, 4, nx, ny, nz,
                                false, false, false);

  int pos = 0, i, j, k;
  bbox_t bb;
  int* this_tuple = int_tuple_new(3);
  int* prev_tuple = NULL;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, &bb))
  {
    this_tuple[0] = i;
    this_tuple[1] = j;
    this_tuple[2] = k;
    if (prev_tuple == NULL)
    {
      prev_tuple = int_tuple_new(3);
      prev_tuple[0] = i;
      prev_tuple[1] = j;
      prev_tuple[2] = k;
    }
    else
    {
      assert_int_equal(-1, int_tuple_cmp(prev_tuple, this_tuple));
      int_tuple_copy(this_tuple, prev_tuple);
    }
  }
  int_tuple_free(this_tuple);
  if (prev_tuple != NULL)
    int_tuple_free(prev_tuple);
  unimesh_free(mesh);
}

static void test_repartition(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};

  // This is created using the naive partitioning.
  unimesh_t* mesh = unimesh_new(MPI_COMM_WORLD, &bbox, 4, 4, 4, nx, ny, nz,
                                false, false, false);

  // Repartition it!
  repartition_unimesh(&mesh, NULL, 0.05, NULL, 0);

  // Just a smoke test for now, folks.
  unimesh_free(mesh);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_ctors),
    cmocka_unit_test(test_next_patch),
    cmocka_unit_test(test_repartition)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
