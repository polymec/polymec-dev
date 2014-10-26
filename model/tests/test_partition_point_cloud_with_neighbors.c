// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
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
  exchanger_t* distributor = partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05);
  exchanger_free(distributor);

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

  // Plot it.
  double p[cloud->num_points];
  for (int i = 0; i < cloud->num_points; ++i)
    p[i] = 1.0*rank;
  char filename[FILENAME_MAX];
  snprintf(filename, FILENAME_MAX, "linear_cloud_partition_with_neighbors_%d", N);
  silo_file_t* silo = silo_file_new(comm, filename, filename, 1, 0, 0, 0.0);
  silo_file_write_point_cloud(silo, "cloud", cloud);
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p);
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
  real_t dx = 1.0/ny, dy = 1.0/ny;
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
  exchanger_t* distributor = partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05);
  exchanger_free(distributor);

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
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p);
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
  real_t dx = 1.0/ny, dy = 1.0/ny, dz = 1.0/nz;
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
  exchanger_t* distributor = partition_point_cloud_with_neighbors(&cloud, &pairing, comm, NULL, 0.05);
  exchanger_free(distributor);

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
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p);
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
  test_partition_cubic_cloud(state, 50, 50, 50);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_partition_small_linear_cloud),
    unit_test(test_partition_large_linear_cloud),
    unit_test(test_partition_small_planar_cloud),
    unit_test(test_partition_large_planar_cloud),
    unit_test(test_partition_small_cubic_cloud),
    unit_test(test_partition_large_cubic_cloud)
  };
  return run_tests(tests);
}
