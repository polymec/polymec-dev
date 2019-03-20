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
#include "geometry/partition_polymesh.h"
#include "geometry/create_uniform_polymesh.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
static const real_t volume_tolerance = 1e-12;
static const real_t area_tolerance = 1e-12;
#else
static const real_t volume_tolerance = 1e-6;
static const real_t area_tolerance = 1e-6;
#endif
static void test_partition_linear_mesh(void** state)
{
  // Create a 100x1x1 uniform mesh.
  int nx = 100, ny = 1, nz = 1;
  real_t dx = 1.0/nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = dx, .z1 = 0.0, .z2 = dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  assert_true(partition_polymesh(&mesh, MPI_COMM_WORLD, NULL, 0.05, NULL, 0));

  // Check the ghost cells.
  int rank, nprocs;
  MPI_Comm_rank(mesh->comm, &rank);
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs > 1)
  {
    exchanger_t* ex = polymesh_exchanger(mesh, POLYMESH_CELL);
    int pos = 0, proc, *indices, num_indices;
    int num_sends = 0, num_receives = 0;
    while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
      num_sends += num_indices;
    pos = 0;
    while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
      num_receives += num_indices;
    assert_true((mesh->num_ghost_cells == 1) || (mesh->num_ghost_cells == 2));
    assert_true(num_sends == mesh->num_ghost_cells);
    assert_true(num_receives == mesh->num_ghost_cells);
    ex = NULL;
  }
  else
    assert_int_equal(0, mesh->num_ghost_cells);

  // Check the geometry of the mesh.
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (!reals_nearly_equal(mesh->cell_volumes[c], dx*dx*dx, volume_tolerance))
    {
      cell_volumes_are_ok = 0;
      break;
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (!reals_nearly_equal(mesh->face_areas[f], dx*dx, area_tolerance))
    {
      face_areas_are_ok = 0;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &face_areas_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Check the resulting exchanger.
  exchanger_is_valid(polymesh_exchanger(mesh, POLYMESH_CELL), NULL);

  // Clean up.
  polymesh_free(mesh);
}

static void test_partition_slab_mesh(void** state)
{
  // Create a 10x10x1 uniform mesh.
  int nx = 10, ny = 10, nz = 1;
  real_t dx = 1.0/nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  assert_true(partition_polymesh(&mesh, MPI_COMM_WORLD, NULL, 0.05, NULL, 0));

  // Check the geometry of the mesh.
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (!reals_nearly_equal(mesh->cell_volumes[c], dx*dx*dx, volume_tolerance))
    {
      cell_volumes_are_ok = 0;
      break;
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (!reals_nearly_equal(mesh->face_areas[f], dx*dx, area_tolerance))
    {
      face_areas_are_ok = 0;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &face_areas_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Clean up.
  polymesh_free(mesh);
}

static void test_partition_box_mesh(void** state)
{
  // Create a 10x10x10 uniform mesh.
  int nx = 10, ny = 10, nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  assert_true(partition_polymesh(&mesh, MPI_COMM_WORLD, NULL, 0.05, NULL, 0));

  // Check the geometry of the mesh.
  real_t dx = 1.0/nx;
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (!reals_nearly_equal(mesh->cell_volumes[c], dx*dx*dx, volume_tolerance))
    {
      cell_volumes_are_ok = 0;
      break;
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (!reals_nearly_equal(mesh->face_areas[f], dx*dx, area_tolerance))
    {
      face_areas_are_ok = 0;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &face_areas_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Clean up.
  polymesh_free(mesh);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_partition_linear_mesh),
    cmocka_unit_test(test_partition_slab_mesh),
    cmocka_unit_test(test_partition_box_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
