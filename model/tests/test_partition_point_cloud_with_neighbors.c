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
#include "geometry/create_point_lattice.h"
#include "model/partition_point_cloud_with_neighbors.h"
//#include "io/silo_file.h"

// This creates a neighbor pairing using a hat function.
extern neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h);

static void test_partition_linear_cloud(void** state)
{
  int N = 10;

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
  assert_true(partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05, NULL, 0));

  // Now check data. Points should all fall on dx tick marks.
  for (int i = 0; i < cloud->num_points; ++i)
  {
    real_t x = cloud->points[i].x;
    real_t y = cloud->points[i].y;
    real_t z = cloud->points[i].z;
    int j = (int)lround(x/dx - 0.5);
    assert_true(reals_nearly_equal(x, (0.5+j)*dx, 1e-6));
    assert_true(reals_nearly_equal(y, 0.5, 1e-6));
    assert_true(reals_nearly_equal(z, 0.5, 1e-6));
  }

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);
}

static void test_partition_planar_cloud(void** state)
{
  int nx = 10, ny = 10;

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
  assert_true(partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05, NULL, 0));

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);
}

static void test_partition_cubic_cloud(void** state)
{
  int nx = 10, ny = 10, nz = 10;

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
  assert_true(partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05, NULL, 0));

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_partition_linear_cloud),
    cmocka_unit_test(test_partition_planar_cloud),
    cmocka_unit_test(test_partition_cubic_cloud)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
