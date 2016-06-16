// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "core/silo_file.h"
#include "core/partition_mesh.h"
#include "geometry/create_uniform_mesh.h"

static void test_migrate_mesh_1_proc(void** state)
{
  // Create a mesh.
  int nx = 10, ny = 10, nz = 10;
  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify_topology(mesh, polymec_error);

  int64_t P[mesh->num_cells];
  memset(P, 0, sizeof(int64_t) * mesh->num_cells);
  migrator_t* m = migrate_mesh(&mesh, MPI_COMM_WORLD, P);
  assert_true(m != NULL);
  migrator_free(m);

  mesh_free(mesh);
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
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify_topology(mesh, polymec_error);

  // We'll swap the cells on the processes.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int other = (rank == 1) ? 0 : 1;
  int64_t P[2] = {other, other};
  migrator_t* m = migrate_mesh(&mesh, MPI_COMM_WORLD, P);
  assert_true(m != NULL);
  migrator_free(m);

  // Check the numbers.
  assert_int_equal(2, mesh->num_cells);
  assert_int_equal(11, mesh->num_faces);
  assert_int_equal(20, mesh->num_edges);
  assert_int_equal(12, mesh->num_nodes);
  assert_true(ABS(mesh->cell_volumes[0] - 1.0) < 1e-12);
  assert_true(ABS(mesh->cell_volumes[1] - 1.0) < 1e-12);

  mesh_free(mesh);
}

static void test_migrate_4x1x1_mesh_3_proc(void** state)
{
  // Create a 4x1x1 mesh across 3 processes.
  int nx = 4, ny = 1, nz = 1;

  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify_topology(mesh, polymec_error);

  // Move one cell to the next process.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int64_t P[mesh->num_cells];
  for (int i = 0; i < mesh->num_cells-1; ++i)
    P[i] = rank;
  int right = (rank + 1) % 3;
  P[mesh->num_cells-1] = right;
  migrator_t* m = migrate_mesh(&mesh, MPI_COMM_WORLD, P);
  assert_true(m != NULL);
  migrator_free(m);

  // Check the numbers.
  if (rank == 1)
  {
    assert_int_equal(2, mesh->num_cells);
    assert_int_equal(11, mesh->num_faces);
    assert_int_equal(20, mesh->num_edges);
    assert_int_equal(12, mesh->num_nodes);
    assert_true(ABS(mesh->cell_volumes[0] - 1.0) < 1e-12);
    assert_true(ABS(mesh->cell_volumes[1] - 1.0) < 1e-12);
  }
  else
  {
    assert_int_equal(1, mesh->num_cells);
    assert_int_equal(6, mesh->num_faces);
    assert_int_equal(12, mesh->num_edges);
    assert_int_equal(8, mesh->num_nodes);
    assert_true(ABS(mesh->cell_volumes[0] - 1.0) < 1e-12);
  }

  mesh_free(mesh);
}

static void test_migrate_4x1x1_mesh_4_proc(void** state)
{
  // Create a 4x1x1 mesh across 4 processes.
  int nx = 4, ny = 1, nz = 1;

  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify_topology(mesh, polymec_error);

  // Cycle our cell to the previous process.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int left = (rank == 0) ? 3 : (rank - 1) % 4;
  int64_t P[1] = {left};
  migrator_t* m = migrate_mesh(&mesh, MPI_COMM_WORLD, P);
  assert_true(m != NULL);
  migrator_free(m);

  // Check the numbers.
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);
  assert_true(ABS(mesh->cell_volumes[0] - 1.0) < 1e-12);

  mesh_free(mesh);
}

TEST_MIGRATE_MESH(test_migrate_4x1x1_mesh)

static void test_repartition_uniform_mesh_of_size(void** state, int nx, int ny, int nz)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create an nx x ny x nz uniform mesh.
  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify_topology(mesh, polymec_error);

  // Repartition it.
  migrator_t* m = repartition_mesh(&mesh, NULL, 0.05);
  migrator_verify(m, polymec_error);

for (int n = 0; n < mesh->num_nodes; ++n)
log_debug("node %d: (%g, %g, %g)", n, mesh->nodes[n].x, mesh->nodes[n].y, mesh->nodes[n].z);

  // Since the mesh is uniform, we can check the properties of each cell.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    assert_int_equal(6, mesh_cell_num_faces(mesh, c));
    real_t V = dx * dx * dx;
log_debug("cell %d at (%g, %g, %g): V = %g, should be %g", c, mesh->cell_centers[c].x, mesh->cell_centers[c].y, mesh->cell_centers[c].z, mesh->cell_volumes[c], V);
    assert_true(ABS(mesh->cell_volumes[c] - V)/V < 1e-14);
  }

  // We can also check the properties of each face.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    assert_int_equal(4, mesh_face_num_edges(mesh, f));
    assert_int_equal(4, mesh_face_num_nodes(mesh, f));
    real_t A = dx * dx;
    assert_true(ABS(mesh->face_areas[f] - A)/A < 1e-14);
  }

  // Check the resulting exchanger.
  exchanger_verify(mesh_exchanger(mesh), polymec_error);

  // Plot it.
  double r[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    r[c] = 1.0*rank;
  char prefix[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "%dx%dx%d_uniform_mesh_repartition", nx, ny, nz);
  silo_file_t* silo = silo_file_new(mesh->comm, prefix, prefix, 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", r, NULL);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_repartition_4x1x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 1, 1);
}

void test_repartition_2x2x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 2, 2, 1);
}

void test_repartition_4x4x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 4, 1);
}

void test_repartition_2x2x2_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 2, 2, 2);
}

void test_repartition_4x4x4_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 4, 4);
}

void test_repartition_128x1x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 128, 1, 1);
}

void test_repartition_128x128x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 128, 128, 1);
}

void test_repartition_32x32x32_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 32, 32, 32);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_migrate_4x1x1_mesh),
    cmocka_unit_test(test_repartition_4x1x1_uniform_mesh),
    cmocka_unit_test(test_repartition_2x2x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_4x4x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_2x2x2_uniform_mesh),
//    cmocka_unit_test(test_repartition_4x4x4_uniform_mesh),
//    cmocka_unit_test(test_repartition_128x1x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_128x128x1_uniform_mesh),
//    cmocka_unit_test(test_repartition_32x32x32_uniform_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
