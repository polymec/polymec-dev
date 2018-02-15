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
#include "geometry/unimesh_patch_bc.h"

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

static void set_up_bcs_if_needed(unimesh_field_t* field)
{
  unimesh_t* mesh = unimesh_field_mesh(field);
  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);

  // If the mesh isn't periodic, we set up zero boundary conditions on 
  // the field at the boundaries of the mesh.
  real_t zeros[3] = {0.0, 0.0, 0.0};
  unimesh_patch_bc_t* zero_bc = constant_unimesh_patch_bc_new(mesh, zeros, 3);

  int pos = 0, pi, pj, pk;
  while (unimesh_next_patch(mesh, &pos, &pi, &pj, &pk, NULL))
  {
    if (!x_periodic)
    {
      if (pi == 0)
        unimesh_field_set_patch_bc(field, pi, pj, pk, UNIMESH_X1_BOUNDARY, zero_bc);
      else if (pi == npx-1)
        unimesh_field_set_patch_bc(field, pi, pj, pk, UNIMESH_X2_BOUNDARY, zero_bc);
    }
    if (!y_periodic)
    {
      if (pj == 0)
        unimesh_field_set_patch_bc(field, pi, pj, pk, UNIMESH_Y1_BOUNDARY, zero_bc);
      else if (pj == npy-1)
        unimesh_field_set_patch_bc(field, pi, pj, pk, UNIMESH_Y2_BOUNDARY, zero_bc);
    }
    if (!z_periodic)
    {
      if (pk == 0)
        unimesh_field_set_patch_bc(field, pi, pj, pk, UNIMESH_Z1_BOUNDARY, zero_bc);
      else if (pk == npz-1)
        unimesh_field_set_patch_bc(field, pi, pj, pk, UNIMESH_Z2_BOUNDARY, zero_bc);
    }
  }
  polymec_release(zero_bc);
}

static void test_cell_field(void** state, unimesh_t* mesh)
{
  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  unimesh_field_t* field = unimesh_field_new(mesh, UNIMESH_CELL, 3);

  set_up_bcs_if_needed(field);

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
    DECLARE_UNIMESH_CELL_ARRAY(f, patch);

    // interior values
    for (int i = 1; i <= patch->nx; ++i)
    {
      for (int j = 1; j <= patch->ny; ++j)
      {
        for (int k = 1; k <= patch->nz; ++k)
        {
           assert_true(reals_equal(f[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(f[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(f[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    int pi_m = (pi > 0) ? pi - 1 
                        : x_periodic ? npx-1 : 0;
    int pi_p = (pi < npx-1) ? pi + 1 : 0;
    for (int j = 1; j <= patch->ny; ++j)
    {
      for (int k = 1; k <= patch->nz; ++k)
      {
        assert_true(reals_equal(f[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(f[patch->nx+1][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries.
    int pj_m = (pj > 0) ? pj - 1 
                        : y_periodic ? npy-1 : 0;
    int pj_p = (pj < npy-1) ? pj + 1 : 0;
    for (int i = 1; i <= patch->nx; ++i)
    {
      for (int k = 1; k <= patch->nz; ++k)
      {
        assert_true(reals_equal(f[i][0][k][1], 1.0 * pj_m));
        assert_true(reals_equal(f[i][patch->ny+1][k][1], 1.0 * pj_p));
      }
    }

    // z boundaries.
    int pk_m = (pk > 0) ? pk - 1 
                        : z_periodic ? npz-1 : 0;
    int pk_p = (pk < npz-1) ? pk + 1 : 0;
    for (int i = 1; i <= patch->nx; ++i)
    {
      for (int j = 1; j <= patch->ny; ++j)
      {
        assert_true(reals_equal(f[i][j][0][2], 1.0 * pk_m));
        assert_true(reals_equal(f[i][j][patch->nz+1][2], 1.0 * pk_p));
      }
    }
  }

  // Clean up.
  unimesh_field_free(field);
  unimesh_free(mesh);
}

static void test_node_field(void** state, unimesh_t* mesh)
{
  int npx, npy, npz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  bool x_periodic, y_periodic, z_periodic;
  unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
  unimesh_field_t* field = unimesh_field_new(mesh, UNIMESH_NODE, 3);

  set_up_bcs_if_needed(field);

  // Fill our field with patch-specific values.
  int pos = 0, pi, pj, pk;
  unimesh_patch_t* patch;
  while (unimesh_field_next_patch(field, &pos, &pi, &pj, &pk, &patch, NULL))
  {
    assert_true(patch->centering == UNIMESH_NODE);
    assert_true(patch->nc == 3);
    DECLARE_UNIMESH_NODE_ARRAY(f, patch);
    for (int i = 0; i <= patch->nx; ++i)
    {
      for (int j = 0; j <= patch->ny; ++j)
      {
        for (int k = 0; k <= patch->nz; ++k)
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
    DECLARE_UNIMESH_NODE_ARRAY(f, patch);

    // interior values
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        for (int k = 1; k < patch->nz; ++k)
        {
           assert_true(reals_equal(f[i][j][k][0], 1.0 * pi));
           assert_true(reals_equal(f[i][j][k][1], 1.0 * pj));
           assert_true(reals_equal(f[i][j][k][2], 1.0 * pk));
        }
      }
    }

    // x boundaries
    int pi_m = (pi > 0) ? pi - 1 : 0;
    int pi_p = (pi < npx-1) ? pi : 0;
    for (int j = 1; j < patch->ny; ++j)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(f[0][j][k][0], 1.0 * pi_m));
        assert_true(reals_equal(f[patch->nx][j][k][0], 1.0 * pi_p));
      }
    }

    // y boundaries 
    int pj_m = (pj > 0) ? pj - 1 : 0;
    int pj_p = (pj < npy-1) ? pj : 0;
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int k = 1; k < patch->nz; ++k)
      {
        assert_true(reals_equal(f[i][0][k][1], 1.0 * pj_m));
        assert_true(reals_equal(f[i][patch->ny][k][1], 1.0 * pj_p));
      }
    }

    // z boundaries 
    int pk_m = (pk > 0) ? pk - 1 : 0;
    int pk_p = (pk < npz-1) ? pk : 0;
    for (int i = 1; i < patch->nx; ++i)
    {
      for (int j = 1; j < patch->ny; ++j)
      {
        assert_true(reals_equal(f[i][j][0][2], 1.0 * pk_m));
        assert_true(reals_equal(f[i][j][patch->nz][2], 1.0 * pk_p));
      }
    }
  }

  // Clean up.
  unimesh_field_free(field);
  unimesh_free(mesh);
}

static void test_serial_periodic_cell_field(void** state)
{
  unimesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

static void test_serial_periodic_face_fields(void** state)
{
}

static void test_serial_periodic_edge_fields(void** state)
{
}

static void test_serial_periodic_node_field(void** state)
{
  unimesh_t* mesh = periodic_mesh(MPI_COMM_SELF);
  test_node_field(state, mesh);
}

static void test_serial_aperiodic_cell_field(void** state)
{
  unimesh_t* mesh = aperiodic_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

static void test_serial_aperiodic_face_fields(void** state)
{
}

static void test_serial_aperiodic_edge_fields(void** state)
{
}

static void test_serial_aperiodic_node_field(void** state)
{
  unimesh_t* mesh = aperiodic_mesh(MPI_COMM_SELF);
  test_node_field(state, mesh);
}

static void test_parallel_periodic_cell_field(void** state)
{
  unimesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

static void test_parallel_periodic_face_fields(void** state)
{
}

static void test_parallel_periodic_edge_fields(void** state)
{
}

static void test_parallel_periodic_node_field(void** state)
{
  unimesh_t* mesh = periodic_mesh(MPI_COMM_WORLD);
  test_node_field(state, mesh);
}

static void test_parallel_aperiodic_cell_field(void** state)
{
  unimesh_t* mesh = aperiodic_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

static void test_parallel_aperiodic_face_fields(void** state)
{
}

static void test_parallel_aperiodic_edge_fields(void** state)
{
}

static void test_parallel_aperiodic_node_field(void** state)
{
  unimesh_t* mesh = aperiodic_mesh(MPI_COMM_WORLD);
  test_node_field(state, mesh);
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
