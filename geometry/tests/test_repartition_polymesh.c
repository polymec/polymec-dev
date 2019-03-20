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
static const real_t abs_tolerance = 1e-12;
static const real_t frac_tolerance = 1e-14;
#else
static const real_t abs_tolerance = 1e-6;
static const real_t frac_tolerance = 1e-6;
#endif

static void test_migrate_mesh_1_proc(void** state)
{
  // Create a mesh.
  int nx = 10, ny = 10, nz = 10;
  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  assert_true(polymesh_is_valid(mesh, NULL));

  int64_t P[mesh->num_cells];
  memset(P, 0, sizeof(int64_t) * mesh->num_cells);
  redistribute_polymesh(&mesh, P, NULL, 0);

  polymesh_free(mesh);
}

#define TEST_MIGRATE_MESH(function_name) \
static void function_name(void** state) \
{ \
  int nprocs; \
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); \
  switch (nprocs) \
  { \
    case 1: \
      test_migrate_mesh_1_proc(state); \
      break; \
    case 2: \
      function_name##_2_proc(state); \
      break; \
    case 3: \
      function_name##_3_proc(state); \
      break; \
    case 4: \
      function_name##_4_proc(state); \
      break; \
  } \
}

static void test_migrate_4x1x1_mesh_2_proc(void** state)
{
  // Create a 4x1x1 mesh across 2 processes.
  int nx = 4, ny = 1, nz = 1;

  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  assert_true(polymesh_is_valid(mesh, NULL));

  // We'll swap the cells on the processes.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int other = (rank == 1) ? 0 : 1;
  int64_t P[2] = {other, other};
  redistribute_polymesh(&mesh, P, NULL, 0);

  // Check the numbers.
  assert_int_equal(2, mesh->num_cells);
  assert_int_equal(11, mesh->num_faces);
  assert_int_equal(20, mesh->num_edges);
  assert_int_equal(12, mesh->num_nodes);
  assert_true(ABS(mesh->cell_volumes[0] - 1.0) < abs_tolerance);
  assert_true(ABS(mesh->cell_volumes[1] - 1.0) < abs_tolerance);

  polymesh_free(mesh);
}

static void test_migrate_4x1x1_mesh_3_proc(void** state)
{
  // Create a 4x1x1 mesh across 3 processes.
  int nx = 4, ny = 1, nz = 1;

  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  assert_true(polymesh_is_valid(mesh, NULL));

  // Move one cell to the next process.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int64_t P[mesh->num_cells];
  for (int i = 0; i < mesh->num_cells-1; ++i)
    P[i] = rank;
  int right = (rank + 1) % 3;
  P[mesh->num_cells-1] = right;
  redistribute_polymesh(&mesh, P, NULL, 0);

  // Check the numbers.
  if (rank == 2)
  {
    assert_int_equal(2, mesh->num_cells);
    assert_int_equal(11, mesh->num_faces);
    assert_int_equal(20, mesh->num_edges);
    assert_int_equal(12, mesh->num_nodes);
    assert_true(ABS(mesh->cell_volumes[0] - 1.0) < abs_tolerance);
    assert_true(ABS(mesh->cell_volumes[1] - 1.0) < abs_tolerance);
  }
  else
  {
    assert_int_equal(1, mesh->num_cells);
    assert_int_equal(6, mesh->num_faces);
    assert_int_equal(12, mesh->num_edges);
    assert_int_equal(8, mesh->num_nodes);
    assert_true(ABS(mesh->cell_volumes[0] - 1.0) < abs_tolerance);
  }

  polymesh_free(mesh);
}

static void test_migrate_4x1x1_mesh_4_proc(void** state)
{
  // Create a 4x1x1 mesh across 4 processes.
  int nx = 4, ny = 1, nz = 1;

  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  assert_true(polymesh_is_valid(mesh, NULL));

  // Cycle our cell to the previous process.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int left = (rank == 0) ? 3 : (rank - 1) % 4;
  int64_t P[1] = {left};
  redistribute_polymesh(&mesh, P, NULL, 0);

  // Check the numbers.
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);
  assert_true(ABS(mesh->cell_volumes[0] - 1.0) < abs_tolerance);

  polymesh_free(mesh);
}

TEST_MIGRATE_MESH(test_migrate_4x1x1_mesh)

static void test_repartition_uniform_mesh_of_size(void** state, int nx, int ny, int nz)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create an nx x ny x nz uniform mesh.
  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  assert_true(polymesh_is_valid(mesh, NULL));

  // Repartition it.
  assert_true(repartition_polymesh(&mesh, NULL, 0.05, NULL, 0));

  // Since the mesh is uniform, we can check the properties of each cell.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    assert_int_equal(6, polymesh_cell_num_faces(mesh, c));
    real_t V = dx * dx * dx;
    assert_true(ABS(mesh->cell_volumes[c] - V)/V < frac_tolerance);
  }

  // We can also check the properties of each face.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    assert_int_equal(4, polymesh_face_num_edges(mesh, f));
    assert_int_equal(4, polymesh_face_num_nodes(mesh, f));
    real_t A = dx * dx;
    assert_true(ABS(mesh->face_areas[f] - A)/A < frac_tolerance);
  }

  // Check the resulting exchanger.
  assert_true(exchanger_is_valid(polymesh_exchanger(mesh, POLYMESH_CELL), NULL));

  // Clean up.
  polymesh_free(mesh);
}

static void test_repartition_4x1x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 1, 1);
}

static void test_repartition_2x2x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 2, 2, 1);
}

static void test_repartition_4x4x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 4, 1);
}

#if 0
static void test_repartition_2x2x2_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 2, 2, 2);
}

static void test_repartition_4x4x4_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 4, 4);
}

static void test_repartition_128x1x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 128, 1, 1);
}

static void test_repartition_128x128x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 128, 128, 1);
}

static void test_repartition_32x32x32_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 32, 32, 32);
}
#endif

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_migrate_4x1x1_mesh),
    cmocka_unit_test(test_repartition_4x1x1_uniform_mesh),
    cmocka_unit_test(test_repartition_2x2x1_uniform_mesh),
    cmocka_unit_test(test_repartition_4x4x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_2x2x2_uniform_mesh),
//    cmocka_unit_test(test_repartition_4x4x4_uniform_mesh),
//    cmocka_unit_test(test_repartition_128x1x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_128x128x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_32x32x32_uniform_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
