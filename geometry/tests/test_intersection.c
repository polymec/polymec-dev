// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "geometry/plane.h"
#include "geometry/intersection.h"
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
  planes[0] = plane_new(&n1, &x1);
  planes[1] = plane_new(&n2, &x2);
  planes[2] = plane_new(&n3, &x3);
  planes[3] = plane_new(&n4, &x4);
  planes[4] = plane_new(&n5, &x5);
  planes[5] = plane_new(&n6, &x6);

  // Now construct their intersection, which should be a cube.
  sp_func_t* i = intersection_new(planes, 6);
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
  planes[0] = plane_new(&n1, &x1);
  planes[1] = plane_new(&n2, &x2);
  planes[2] = plane_new(&n3, &x3);
  planes[3] = plane_new(&n4, &x4);
  planes[4] = plane_new(&n5, &x5);
  planes[5] = plane_new(&n6, &x6);
  sp_func_t* i = intersection_new(planes, 6);
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  generate_octave_script_for_surface(i, 40, &bbox, "test_intersection.m");
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
