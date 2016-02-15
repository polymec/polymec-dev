// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "geometry/plane_sp_func.h"
#include "geometry/intersection_sp_func.h"
#include "generate_octave_script_for_surface.h"

void test_construct(void** state)
{
  // Create six planes.
  vector_t n1 = { 1.0, 0.0, 0.0}, n2 = {-1.0, 0.0, 0.0},
           n3 = { 0.0, 1.0, 0.0}, n4 = { 0.0,-1.0, 0.0},
           n5 = { 0.0, 0.0, 1.0}, n6 = { 0.0, 0.0,-1.0};
  point_t x1 = {-0.5, 0.0, 0.0}, x2 = { 0.5, 0.0, 0.0},
          x3 = { 0.0,-0.5, 0.0}, x4 = { 0.0, 0.5, 0.0},
          x5 = { 0.0, 0.0,-0.5}, x6 = { 0.0, 0.0, 0.5};
  sp_func_t* planes[6];
  planes[0] = plane_sp_func_new(&n1, &x1);
  planes[1] = plane_sp_func_new(&n2, &x2);
  planes[2] = plane_sp_func_new(&n3, &x3);
  planes[3] = plane_sp_func_new(&n4, &x4);
  planes[4] = plane_sp_func_new(&n5, &x5);
  planes[5] = plane_sp_func_new(&n6, &x6);

  // Now construct their intersection, which should be a cube.
  sp_func_t* i = intersection_sp_func_new(planes, 6);
  assert_true(sp_func_num_comp(i) == 1);
  assert_true(sp_func_has_deriv(i, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  vector_t n1 = { 1.0, 0.0, 0.0}, n2 = {-1.0, 0.0, 0.0},
           n3 = { 0.0, 1.0, 0.0}, n4 = { 0.0,-1.0, 0.0},
           n5 = { 0.0, 0.0, 1.0}, n6 = { 0.0, 0.0,-1.0};
  point_t x1 = {-0.5, 0.0, 0.0}, x2 = { 0.5, 0.0, 0.0},
          x3 = { 0.0,-0.5, 0.0}, x4 = { 0.0, 0.5, 0.0},
          x5 = { 0.0, 0.0,-0.5}, x6 = { 0.0, 0.0, 0.5};
  sp_func_t* planes[6];
  planes[0] = plane_sp_func_new(&n1, &x1);
  planes[1] = plane_sp_func_new(&n2, &x2);
  planes[2] = plane_sp_func_new(&n3, &x3);
  planes[3] = plane_sp_func_new(&n4, &x4);
  planes[4] = plane_sp_func_new(&n5, &x5);
  planes[5] = plane_sp_func_new(&n6, &x6);
  sp_func_t* i = intersection_sp_func_new(planes, 6);
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  generate_octave_script_for_surface(i, 40, &bbox, "test_intersection.m");
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_construct),
    cmocka_unit_test(test_plot)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
