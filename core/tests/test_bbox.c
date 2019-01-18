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
#include "core/polymec.h"
#include "core/point.h"

static void test_ctor(void** state)
{
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  assert_true(reals_equal(bbox->x1, 0.0));
  assert_true(reals_equal(bbox->x2, 1.0));
  assert_true(reals_equal(bbox->y1, 0.0));
  assert_true(reals_equal(bbox->y2, 1.0));
  assert_true(reals_equal(bbox->z1, 0.0));
  assert_true(reals_equal(bbox->z2, 1.0));
  bbox_fprintf(bbox, stdout);

  assert_false(bbox_is_empty_set(bbox));
  assert_false(bbox_is_point(bbox));
  assert_false(bbox_is_line(bbox));
  assert_false(bbox_is_plane(bbox));
  bbox_t* bbox2 = bbox_clone(bbox);
  assert_true(bbox_intersects_bbox(bbox, bbox2));

  bbox_t* empty = empty_set_bbox_new();
  assert_true(bbox_is_empty_set(empty));
  assert_false(bbox_intersects_bbox(bbox, empty));
}

static void test_intersection(void** state)
{
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  bbox_t* bbox2 = bbox_new(0.5, 1.0, 0.5, 1.0, 0.5, 1.0);
  bbox_t intersection;
  bbox_intersect_bbox(bbox, bbox2, &intersection);
  assert_true(reals_equal(intersection.x1, 0.5));
  assert_true(reals_equal(intersection.x2, 1.0));
  assert_true(reals_equal(intersection.y1, 0.5));
  assert_true(reals_equal(intersection.y2, 1.0));
  assert_true(reals_equal(intersection.z1, 0.5));
  assert_true(reals_equal(intersection.z2, 1.0));
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_intersection)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
