// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
#include "core/octree.h"

static void test_insert(void** state) 
{ 
  // Create a tree containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  octree_t* tree = octree_new(&bounding_box);
  for (int i = 0; i < N; ++i)
  {
    point_t x;
    point_randomize(&x, rng, &bounding_box);
    octree_insert(tree, &x, i);
  }

  assert_int_equal(100, octree_size(tree));
  octree_free(tree);
}

static void test_delete(void** state) 
{ 
  // Create a tree containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  octree_t* tree = octree_new(&bounding_box);
  point_t points[N];
  for (int i = 0; i < N; ++i)
  {
    point_randomize(&points[i], rng, &bounding_box);
    octree_insert(tree, &points[i], i);
  }

  // Now delete each of them individually.
  for (int i = 0; i < N; ++i) 
  {
    octree_delete(tree, i);
    assert_int_equal(N-i-1, octree_size(tree));
  }

  octree_free(tree);
}

static void test_find_nearest(void** state) 
{ 
  // Create a tree containing 100 random points.
  rng_t* rng = host_rng_new();
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  octree_t* tree = octree_new(&bounding_box);
  point_t points[N];
  for (int i = 0; i < N; ++i)
  {
    point_randomize(&points[i], rng, &bounding_box);
    octree_insert(tree, &points[i], i);
  }

  // Now do some nearest point queries and check the answers.
  for (int i = 0; i < 10; ++i) // 10 queries.
  {
    // Pick a random point p.
    point_t p;
    point_randomize(&p, rng, &bounding_box);

    // Which point in the set is p closest to?
    int j = octree_nearest(tree, &p); // Point set says: "j."

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

  octree_free(tree);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_insert),
    cmocka_unit_test(test_delete),
    cmocka_unit_test(test_find_nearest)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
