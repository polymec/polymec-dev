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
#include "core/rng.h"
#include "core/array_utils.h"
#include "geometry/partition_point_cloud.h"

static void test_repartition_linear_cloud(void** state, 
                                          bool balanced, 
                                          bool random_points,
                                          bool weighted)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  real_t imbalance_tol = 0.05;

  // Create a 100x1x1 uniform point cloud distributed over N processes.
  int Np;
  if (nprocs > 1)
  {
    if (balanced)
      Np = 100;
    else
    {
      int N_even = 100, N_odd = 50;
      if ((rank % 2) != 0)
        Np = N_odd;
      else
        Np = N_even;
    }
  }
  else
    Np = 100;

  int N;
  MPI_Allreduce(&Np, &N, 1, MPI_INT, MPI_SUM, comm);
  real_t dx = 1.0/N;
  point_cloud_t* cloud = point_cloud_new(comm, Np);
  if (random_points)
  {
    rng_t* rng = host_rng_new();
    for (int i = 0; i < Np; ++i)
      cloud->points[i].x = rng_uniform(rng);
  }
  else
  {
    for (int i = 0; i < Np; ++i)
      cloud->points[i].x = (0.5+i)*dx + (i % nprocs) * 1.0/nprocs;
  }

  // Repartition it.
  if (!weighted)
  {
    assert_true(repartition_point_cloud(&cloud, NULL, imbalance_tol, NULL, 0));

    // Check the number of points on each domain. 
    assert_true(ABS(1.0*(nprocs*cloud->num_points - N)/N) < imbalance_tol);
  }
  else
  {
    // Provide an imbalance to work with.
    int loads[Np];
    point_cloud_field_t* f_loads = point_cloud_field_new(cloud, 1);
    if (rank % 2)
    {
      int_fill(loads, Np, 10);
      real_fill(f_loads->data, Np, 10.0);
    }
    else
    {
      int_fill(loads, Np, 12);
      real_fill(f_loads->data, Np, 12.0);
    }
    assert_true(repartition_point_cloud(&cloud, loads, imbalance_tol, &f_loads, 1));

    // Check that the loads are balanced.
    assert_int_equal(cloud->num_points, f_loads->num_local_values);
    real_t my_load = 0.0;
    for (int i = 0; i < cloud->num_points; ++i)
      my_load += f_loads->data[i];
    real_t total_load;
    MPI_Allreduce(&my_load, &total_load, 1, MPI_REAL_T, MPI_SUM, comm);
    real_t ideal_load = total_load / nprocs;
    assert_true(ABS((my_load - ideal_load))/ideal_load < imbalance_tol);
  }

  // Now check data. 
  if (!random_points) // Points should all fall on dx tick marks.
  {
    for (int i = 0; i < cloud->num_points; ++i)
    {
      real_t x = cloud->points[i].x;
      real_t y = cloud->points[i].y;
      real_t z = cloud->points[i].z;
      int j = (int)lround(x/dx - 0.5);
      if (j < cloud->num_points)
      {
        // j and i are the same if the point originated on this process, 
        // different if it was transferred by the repartitioning. Either
        // way, the distance between them should be an integer multiple
        // of dx.
        real_t L = ABS(x - (0.5+j)*dx);
        assert_true((round(L/dx) - L/dx) < 1e-6);
        assert_true(ABS(y) < 1e-6);
        assert_true(ABS(z) < 1e-6);
      }
    }
  }
  else // Points should all still fall within [0,1] on the x axis.
  {
    for (int i = 1; i < cloud->num_points; ++i)
    {
      real_t x = cloud->points[i].x;
      real_t y = cloud->points[i].y;
      real_t z = cloud->points[i].z;
      assert_true((x >= 0.0) && (x <= 1.0));
      assert_true(ABS(y) < 1e-6);
      assert_true(ABS(z) < 1e-6);
    }
  }

  // Plot it.
  real_t p[cloud->num_points];
  for (int i = 0; i < cloud->num_points; ++i)
    p[i] = 1.0*rank;
  char balance_str[33];
  if (balanced)
    sprintf(balance_str, "balanced");
  else
    sprintf(balance_str, "unbalanced");
  char random_str[33];
  if (random_points)
    sprintf(random_str, "random");
  else
    sprintf(random_str, "uniform");
  char weight_str[33];
  if (weighted)
    sprintf(weight_str, "weighted");
  else
    sprintf(weight_str, "unweighted");
  char dataset_name[1025];
  snprintf(dataset_name, 1024, "%s_%s_%s_linear_cloud_repartition", 
           balance_str, random_str, weight_str);

  // Clean up.
  point_cloud_free(cloud);
}

static void test_repartition_balanced_uniform_unweighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, true, false, false);
}

static void test_repartition_unbalanced_uniform_unweighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, false, false, false);
}

static void test_repartition_balanced_random_unweighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, true, true, false);
}

static void test_repartition_unbalanced_random_unweighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, false, true, false);
}

static void test_repartition_balanced_uniform_weighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, true, false, true);
}

static void test_repartition_unbalanced_uniform_weighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, false, false, true);
}

static void test_repartition_balanced_random_weighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, true, true, true);
}

static void test_repartition_unbalanced_random_weighted_linear_cloud(void** state)
{
  test_repartition_linear_cloud(state, false, true, true);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_repartition_balanced_uniform_unweighted_linear_cloud),
    cmocka_unit_test(test_repartition_balanced_uniform_weighted_linear_cloud),
    cmocka_unit_test(test_repartition_balanced_random_unweighted_linear_cloud),
    cmocka_unit_test(test_repartition_balanced_random_weighted_linear_cloud),
    cmocka_unit_test(test_repartition_unbalanced_uniform_unweighted_linear_cloud),
    cmocka_unit_test(test_repartition_unbalanced_uniform_weighted_linear_cloud),
    cmocka_unit_test(test_repartition_unbalanced_random_unweighted_linear_cloud),
    cmocka_unit_test(test_repartition_unbalanced_random_weighted_linear_cloud)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
