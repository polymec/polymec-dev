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
#include "geometry/create_point_lattice.h"
#include "model/partition_point_cloud_with_neighbors.h"

// This creates a neighbor pairing using a hat function.
extern neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h);

void test_partition_linear_cloud(void** state, int N)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Create an Nx1x1 uniform point cloud and a neighbor pairing for it.
  real_t dx = 1.0/N;
  point_cloud_t* cloud = NULL;
  neighbor_pairing_t* pairing = NULL;
  if (rank == 0)
  {
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    cloud = create_uniform_point_lattice(MPI_COMM_SELF, N, 1, 1, &bbox);

    // Create a neighbor pairing for the cloud.
    pairing = create_simple_pairing(cloud, 1.2*dx);
  }

  // Partition it.
  migrator_t* distributor = partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05);
  distributor = NULL;

  // Now check data. Points should all fall on dx tick marks.
  for (int i = 0; i < cloud->num_points; ++i)
  {
    real_t x = cloud->points[i].x;
    real_t y = cloud->points[i].y;
    real_t z = cloud->points[i].z;
    int j = lround(x/dx - 0.5);
    assert_true(fabs(x - (0.5+j)*dx) < 1e-6);
    assert_true(fabs(y - 0.5) < 1e-6);
    assert_true(fabs(z - 0.5) < 1e-6);
  }

  // Plot it.
  double p[cloud->num_points];
  for (int i = 0; i < cloud->num_points; ++i)
    p[i] = 1.0*rank;
  char filename[FILENAME_MAX];
  snprintf(filename, FILENAME_MAX, "linear_cloud_partition_with_neighbors_%d", N);
  silo_file_t* silo = silo_file_new(comm, filename, filename, 1, 0, 0, 0.0);
  silo_file_write_point_cloud(silo, "cloud", cloud);
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p, NULL);
  silo_file_close(silo);

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query(filename, filename, &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
}

void test_partition_planar_cloud(void** state, int nx, int ny)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Create a nx x ny x 1 uniform point cloud and a neighbor pairing for it.
  real_t dx = 1.0/nx;
  point_cloud_t* cloud = NULL;
  neighbor_pairing_t* pairing = NULL;
  if (rank == 0)
  {
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    cloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, 1, &bbox);

    // Create a neighbor pairing for the cloud.
    pairing = create_simple_pairing(cloud, 1.2*dx);
  }

  // Partition it.
  migrator_t* distributor = partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05);
  distributor = NULL;

#if 0
  // Now check data. Points should all fall on dx tick marks.
  for (int i = 0; i < cloud->num_points; ++i)
  {
    real_t x = cloud->points[i].x;
    real_t y = cloud->points[i].y;
    real_t z = cloud->points[i].z;
    int j = lround(x/dx - 0.5);
    assert_true(fabs(x - (0.5+j)*dx) < 1e-6);
    assert_true(fabs(y) < 1e-6);
    assert_true(fabs(z) < 1e-6);
  }
#endif

  // Plot it.
  double p[cloud->num_points];
  for (int i = 0; i < cloud->num_points; ++i)
    p[i] = 1.0*rank;
  char filename[FILENAME_MAX];
  snprintf(filename, FILENAME_MAX, "planar_cloud_partition_with_neighbors_%dx%d", nx, ny);
  silo_file_t* silo = silo_file_new(comm, filename, filename, 1, 0, 0, 0.0);
  silo_file_write_point_cloud(silo, "cloud", cloud);
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p, NULL);
  silo_file_close(silo);

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query(filename, filename, &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
}

void test_partition_cubic_cloud(void** state, int nx, int ny, int nz)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Create a nx x ny x 1 uniform point cloud and a neighbor pairing for it.
  real_t dx = 1.0/nx;
  point_cloud_t* cloud = NULL;
  neighbor_pairing_t* pairing = NULL;
  if (rank == 0)
  {
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    cloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, nz, &bbox);

    // Create a neighbor pairing for the cloud.
    pairing = create_simple_pairing(cloud, 1.2*dx);
  }

  // Partition it.
  migrator_t* distributor = partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05);
  distributor = NULL;

#if 0
  // Now check data. Points should all fall on dx tick marks.
  for (int i = 0; i < cloud->num_points; ++i)
  {
    real_t x = cloud->points[i].x;
    real_t y = cloud->points[i].y;
    real_t z = cloud->points[i].z;
    int j = lround(x/dx - 0.5);
    assert_true(fabs(x - (0.5+j)*dx) < 1e-6);
    assert_true(fabs(y) < 1e-6);
    assert_true(fabs(z) < 1e-6);
  }
#endif

  // Plot it.
  double p[cloud->num_points];
  for (int i = 0; i < cloud->num_points; ++i)
    p[i] = 1.0*rank;
  char filename[FILENAME_MAX];
  snprintf(filename, FILENAME_MAX, "cubic_cloud_partition_with_neighbors_%dx%dx%d", nx, ny, nz);
  silo_file_t* silo = silo_file_new(comm, filename, filename, 1, 0, 0, 0.0);
  silo_file_write_point_cloud(silo, "cloud", cloud);
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p, NULL);
  silo_file_close(silo);

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query(filename, filename, &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
}

void test_partition_small_linear_cloud(void** state)
{
  test_partition_linear_cloud(state, 10);
}

void test_partition_large_linear_cloud(void** state)
{
  test_partition_linear_cloud(state, 1000);
}

void test_partition_small_planar_cloud(void** state)
{
  test_partition_planar_cloud(state, 10, 10);
}

void test_partition_large_planar_cloud(void** state)
{
  test_partition_planar_cloud(state, 200, 200);
}

void test_partition_small_cubic_cloud(void** state)
{
  test_partition_cubic_cloud(state, 10, 10, 10);
}

void test_partition_large_cubic_cloud(void** state)
{
  test_partition_cubic_cloud(state, 20, 20, 20);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_partition_small_linear_cloud),
    cmocka_unit_test(test_partition_large_linear_cloud),
    cmocka_unit_test(test_partition_small_planar_cloud),
    cmocka_unit_test(test_partition_large_planar_cloud),
    cmocka_unit_test(test_partition_small_cubic_cloud),
    cmocka_unit_test(test_partition_large_cubic_cloud)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
