// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
#include "core/partition_point_cloud.h"

void test_repartition_linear_cloud(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // Create a 100x1x1 uniform point cloud distributed over N processes.
  int N = 100, Np;
  if (nprocs > 1) 
    Np = (rank < nprocs-1) ? N/nprocs : N - (rank-1)*N/nprocs;
  else
    Np = N;

  real_t dx = 1.0/N;
  point_cloud_t* cloud = point_cloud_new(comm, N);
  for (int i = 0; i < Np; ++i)
    cloud->points[i].x = (0.5+i)*dx + (i % nprocs) * 1.0/nprocs;

  // Repartition it.
  exchanger_t* migrator = repartition_point_cloud(&cloud, NULL, 0.05);
  exchanger_free(migrator);

  // Check the number of points on each domain.
printf("%d: %d %d %g\n", rank, cloud->num_points, Np, fabs(1.0*(cloud->num_points - Np)/Np));
  assert_true(fabs(1.0*(cloud->num_points - Np)/Np) < 0.05);

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
  silo_file_t* silo = silo_file_new(comm, "linear_cloud_repartition", "linear_cloud_repartition", 1, 0, 0, 0.0);
  silo_file_write_point_cloud(silo, "cloud", cloud);
  silo_file_write_scalar_point_field(silo, "rank", "cloud", p);
  silo_file_close(silo);

  // Clean up.
  point_cloud_free(cloud);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query("linear_cloud_repartition", "linear_cloud_repartition",
                              &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_repartition_linear_cloud)
  };
  return run_tests(tests);
}
