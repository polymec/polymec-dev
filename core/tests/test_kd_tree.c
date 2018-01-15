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
#include "core/array_utils.h"
#include "core/tensor2.h"

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

static void test_within_radius(void** state) 
{ 
  // Create a point set containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rng, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  // Now try out our radius query and measure it against a brute force 
  // calculation.
  real_t radius = 0.25;
  for (int i = 0; i < 10; ++i) // 10 queries.
  {
    // Pick a random point p.
    point_t p;
    point_randomize(&p, rng, &bounding_box);

    // Which points are within our radius?
    int_array_t* candidates = kd_tree_within_radius(tree, &p, radius);

    // Do a linear search for points within the radius.
    int_array_t* actual_points = int_array_new();
    for (int k = 0; k < N; ++k)
    {
      real_t dist = point_distance(&p, &points[k]);
      if (dist < radius)
        int_array_append(actual_points, k);
    }

    // The two arrays must contain the same indices.
    assert_int_equal(actual_points->size, candidates->size);
    for (int k = 0; k < actual_points->size; ++k)
      assert_true(int_lsearch(candidates->data, candidates->size, actual_points->data[k]) != NULL);

    // Clean up.
    int_array_free(candidates);
    int_array_free(actual_points);
  }

  kd_tree_free(tree);
}

// Returns true if y is in the given ellipse centered at x.
static bool point_in_ellipse(void* context, point_t* x, point_t* y)
{
  sym_tensor2_t* E = context;
  vector_t d;
  point_displacement(x, y, &d);
  return (sym_tensor2_ddot(E, &d, &d) < 1.0);
}

static void test_within_ellipse(void** state) 
{ 
  // Create a point set containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rng, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  // Now use a predicate query and compare it to the results of a brute force 
  // calculation.
  sym_tensor2_t E = {0.2, 0.0, 0.0, 
                          0.3, 0.0,
                               0.4};
  for (int i = 0; i < 10; ++i) // 10 queries.
  {
    // Pick a random point p.
    point_t p;
    point_randomize(&p, rng, &bounding_box);

    // Which points fall within the ellipse E at this point?
    int_array_t* candidates = kd_tree_for_predicate(tree, &p, point_in_ellipse, &E);

    // Do a linear search for points within the radius.
    int_array_t* actual_points = int_array_new();
    for (int k = 0; k < N; ++k)
    {
      if (point_in_ellipse(&E, &p, &points[k]))
        int_array_append(actual_points, k);
    }

    // The two arrays must contain the same indices.
    assert_int_equal(actual_points->size, candidates->size);
    for (int k = 0; k < actual_points->size; ++k)
      assert_true(int_lsearch(candidates->data, candidates->size, actual_points->data[k]) != NULL);

    // Clean up.
    int_array_free(candidates);
    int_array_free(actual_points);
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
    cmocka_unit_test(test_within_radius),
    cmocka_unit_test(test_within_ellipse),
    cmocka_unit_test(test_find_ghost_points)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
