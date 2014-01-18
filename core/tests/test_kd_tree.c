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
#include <time.h>
#include "cmockery.h"
#include "core/kd_tree.h"

void test_construct(void** state) 
{ 
  // Create a tree containing 100 random points.
  srand((unsigned)time(NULL));
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rand, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  assert_int_equal(100, kd_tree_size(tree));
  kd_tree_free(tree);
}

void test_find_nearest(void** state) 
{ 
  // Create a point set containing 100 random points.
  srand((unsigned)time(NULL));
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_t points[N];
  for (int i = 0; i < N; ++i)
    point_randomize(&points[i], rand, &bounding_box);
  kd_tree_t* tree = kd_tree_new(points, N);

  // Now do some nearest point queries and check the answers.
  for (int i = 0; i < 10; ++i) // 10 queries.
  {
    // Pick a random point p.
    point_t p;
    point_randomize(&p, rand, &bounding_box);

    // Which point in the set is p closest to?
    int j = kd_tree_nearest(tree, &p); // Point set says: "j."

    // Do a linear search for the closest point.
    double min_dist = FLT_MAX;
    int kk = -1;
    for (int k = 0; k < N; ++k)
    {
      double dist = point_distance(&p, &points[k]);
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

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_construct),
    unit_test(test_find_nearest)
  };
  return run_tests(tests);
}
