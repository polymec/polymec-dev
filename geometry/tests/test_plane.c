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
#include "geometry/plane.h"
#include "generate_octave_script_for_surface.h"

void test_construct(void** state)
{
  // Create a plane.
  vector_t n = {1.0, 1.0, 1.0};
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* p = plane_new(&n, &origin);
  assert_true(sp_func_num_comp(p) == 1);
  assert_true(sp_func_has_deriv(p, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  point_t origin = {0.0, 0.0, 0.0};
  vector_t n = {1.0, 1.0, 1.0};
  sp_func_t* p = plane_new(&n, &origin);
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  generate_octave_script_for_surface(p, 20, &bbox, "test_plane.m");
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_construct),
    unit_test(test_plot)
  };
  return run_tests(tests);
}
