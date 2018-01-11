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
#include <time.h>
#include "cmocka.h"
#include "core/kd_tree.h"

static void test_construct(void** state) 
{ 
  // Create a tree containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rng, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  assert_int_equal(100, kd_tree_size(tree));
  kd_tree_free(tree);
}

static void test_find_nearest(void** state) 
{ 
  // Create a point set containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rng, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  // Now do some nearest point queries and check the answers.
  for (int i = 0; i < 10; ++i) // 10 queries.
  {
    // Pick a random point p.
    point_t p;
    point_randomize(&p, rng, &bounding_box);

    // Which point in the set is p closest to?
    int j = kd_tree_nearest(tree, &p); // Point set says: "j."

    // Do a linear search for the closest point.
    real_t min_dist = REAL_MAX;
    int kk = -1;
    for (int k = 0; k < N; ++k)
    {
      real_t dist = point_distance(&p, &points[k]);
      if (dist < min_dist)
      {
        kk = k;
        min_dist = dist;
      }
    }

    // The two queries must agree with each other.
    assert_true(j == kk);
  }

  kd_tree_free(tree);
}

static void test_find_ghost_points(void** state) 
{
  // Create a point set containing 100 random points within a bounding box 
  // [0, 1] x [0, 1] x [0, 1] on each process.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rng, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  // Now find the points within this same bounding box on other processes.
  // HINT: All points on other processes should be found! This is as much a 
  // HINT: stress test as it is a necessary-but-insufficient test for 
  // HINT: correctness.
  real_t R_max = 0.5;
  exchanger_t* ex = kd_tree_find_ghost_points(tree, MPI_COMM_WORLD, R_max);

  // Make sure we've added all the ghost points to the tree.
  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert_int_equal(kd_tree_size(tree), nprocs * N);

  // Make sure the exchanger has entries in it for each of the other processes.
  assert_true(exchanger_num_sends(ex) == (nprocs-1));
  assert_true(exchanger_num_receives(ex) == (nprocs-1));

  int pos = 0, proc, *indices, num_indices;
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    assert_true(proc != rank);
    assert_int_equal(num_indices, N);
    for (int i = 0; i < num_indices; ++i)
    {
      assert_true(indices[i] < N);
    }
  }
  pos = 0;
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
  {
    assert_true(proc != rank);
    assert_int_equal(num_indices, N);
    for (int i = 0; i < num_indices; ++i)
    {
      assert_true(indices[i] >= N);
    }
  }

  // Clean up.
  ex = NULL;
  kd_tree_free(tree);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_construct),
    cmocka_unit_test(test_find_nearest),
    cmocka_unit_test(test_find_ghost_points)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
