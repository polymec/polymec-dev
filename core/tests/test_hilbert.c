// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/hilbert.h"

void test_ctor(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  hilbert_t* h = hilbert_new(&bbox);
  assert_true(h != NULL);
}

void test_index(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  hilbert_t* h = hilbert_new(&bbox);
  point_t x = {0.0, 0.0, 0.0};
  index_t index = hilbert_index(h, &x);
  assert_true(index == 0);

  x.x = 1.0;
  index = hilbert_index(h, &x);
  assert_true(index == 281474976710655);

  x.y = 1.0;
  index = hilbert_index(h, &x);
  assert_true(index == 150790166094994);

  x.z = 1.0;
  index = hilbert_index(h, &x);
  assert_true(index == 201053554793325);
}

void test_reproduce_point(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  hilbert_t* h = hilbert_new(&bbox);
  point_t x = {1.0, 1.0, 0.0};
  index_t index = hilbert_index(h, &x);
  point_t y;
  hilbert_create_point(h, index, &y);
  assert_true(point_distance(&x, &y) < 1e-12);
}

void test_sort_points(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  hilbert_t* h = hilbert_new(&bbox);
  point_t points[4] = {{0.0, 0.0, 0.0},
                       {1.0, 0.0, 0.0},
                       {1.0, 1.0, 0.0},
                       {1.0, 1.0, 1.0}};
  hilbert_sort_points(h, points, 4);
  assert_true(points[0].x == 0.0);
  assert_true(points[0].y == 0.0);
  assert_true(points[0].z == 0.0);
  assert_true(points[1].x == 1.0);
  assert_true(points[1].y == 1.0);
  assert_true(points[1].z == 0.0);
  assert_true(points[2].x == 1.0);
  assert_true(points[2].y == 1.0);
  assert_true(points[2].z == 1.0);
  assert_true(points[3].x == 1.0);
  assert_true(points[3].y == 0.0);
  assert_true(points[3].z == 0.0);
}


int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_index),
    unit_test(test_reproduce_point),
    unit_test(test_sort_points)
  };
  return run_tests(tests);
}
