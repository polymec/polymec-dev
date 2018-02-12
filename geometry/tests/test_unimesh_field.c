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
#include "geometry/unimesh_field.h"

static unimesh_t* periodic_mesh(MPI_Comm comm)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};
  return unimesh_new(comm, &bbox, 
                     4, 4, 4, 10, 10, 10,
                     true, true, true);
}

static unimesh_t* aperiodic_mesh(MPI_Comm comm)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};
  return unimesh_new(comm, &bbox, 
                     4, 4, 4, 10, 10, 10,
                     false, false, false);
}

static void test_serial_periodic_cell_field(void** state)
{
  unimesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  unimesh_field_t* field = unimesh_field_new(mesh, UNIMESH_CELL, 3);

  // Fill our field with patch-specific values.
  int pos = 0, pi, pj, pk;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    assert_true(patch->centering == UNIMESH_CELL);
    assert_true(patch->nc == 3);
    DECLARE_UNIMESH_CELL_ARRAY(f, patch);
    for (int i = 1; i <= patch->nx; ++i)
    {
      for (int j = 1; j <= patch->ny; ++j)
      {
        for (int k = 1; k <= patch->nz; ++k)
        {
          f[i][j][k][0] = 1.0 * pi;
          f[i][j][k][1] = 1.0 * pj;
          f[i][j][k][2] = 1.0 * pk;
        }
      }
    }
  }

  // Update the patch boundaries in the field.
  unimesh_field_update_patch_boundaries(field, 0.0);

  // Check the patch boundary data.
  pos = 0;
  while (unimesh_field_next_patch(field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    assert_true(patch->centering == UNIMESH_CELL);
    assert_true(patch->nc == 3);
    DECLARE_UNIMESH_CELL_ARRAY(f, patch);

    // x boundaries.
    int pi_m = (pi >= 0) ? pi - 1 : npx-1;
    int pi_p = (pi < npx-1) ? pi + 1 : 0;
    for (int j = 1; j <= patch->ny; ++j)
    {
      for (int k = 1; k <= patch->nz; ++k)
      {
        assert_true(reals_equal(f[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(f[patch->nx][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries.
    int pj_m = (pj >= 0) ? pj - 1 : npy-1;
    int pj_p = (pj < npy-1) ? pj + 1 : 0;
    for (int i = 1; i <= patch->nx; ++i)
    {
      for (int k = 1; k <= patch->nz; ++k)
      {
        for (int c = 0; c <= patch->nc; ++c)
        {
          assert_true(reals_equal(f[i][0][k][1], 1.0 * pj_m));
          assert_true(reals_equal(f[i][patch->ny][k][1], 1.0 * pj_p));
        }
      }
    }

    // z boundaries.
    int pk_m = (pk >= 0) ? pk - 1 : npz-1;
    int pk_p = (pk < npz-1) ? pk + 1 : 0;
    for (int i = 1; i <= patch->nx; ++i)
    {
      for (int j = 1; j <= patch->ny; ++j)
      {
        for (int c = 0; c <= patch->nc; ++c)
        {
          assert_true(reals_equal(f[i][j][0][2], 1.0 * pk_m));
          assert_true(reals_equal(f[i][j][patch->nz][2], 1.0 * pk_p));
        }
      }
    }
  }

  unimesh_field_free(field);
  unimesh_free(mesh);
}

static void test_serial_periodic_face_fields(void** state)
{
}

static void test_serial_periodic_edge_fields(void** state)
{
}

static void test_serial_periodic_node_field(void** state)
{
}

static void test_serial_aperiodic_cell_field(void** state)
{
  unimesh_t* mesh = aperiodic_mesh(MPI_COMM_SELF);
  unimesh_free(mesh);
}

static void test_serial_aperiodic_face_fields(void** state)
{
}

static void test_serial_aperiodic_edge_fields(void** state)
{
}

static void test_serial_aperiodic_node_field(void** state)
{
}

static void test_parallel_periodic_cell_field(void** state)
{
  unimesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  unimesh_free(mesh);
}

static void test_parallel_periodic_face_fields(void** state)
{
}

static void test_parallel_periodic_edge_fields(void** state)
{
}

static void test_parallel_periodic_node_field(void** state)
{
}

static void test_parallel_aperiodic_cell_field(void** state)
{
  unimesh_t* mesh = aperiodic_mesh(MPI_COMM_WORLD);
  unimesh_free(mesh);
}

static void test_parallel_aperiodic_face_fields(void** state)
{
}

static void test_parallel_aperiodic_edge_fields(void** state)
{
}

static void test_parallel_aperiodic_node_field(void** state)
{
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_periodic_cell_field),
    cmocka_unit_test(test_serial_periodic_face_fields),
    cmocka_unit_test(test_serial_periodic_edge_fields),
    cmocka_unit_test(test_serial_periodic_node_field),
    cmocka_unit_test(test_serial_aperiodic_cell_field),
    cmocka_unit_test(test_serial_aperiodic_face_fields),
    cmocka_unit_test(test_serial_aperiodic_edge_fields),
    cmocka_unit_test(test_serial_aperiodic_node_field),
    cmocka_unit_test(test_parallel_periodic_cell_field),
    cmocka_unit_test(test_parallel_periodic_face_fields),
    cmocka_unit_test(test_parallel_periodic_edge_fields),
    cmocka_unit_test(test_parallel_periodic_node_field),
    cmocka_unit_test(test_parallel_aperiodic_cell_field),
    cmocka_unit_test(test_parallel_aperiodic_face_fields),
    cmocka_unit_test(test_parallel_aperiodic_edge_fields),
    cmocka_unit_test(test_parallel_aperiodic_node_field),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
