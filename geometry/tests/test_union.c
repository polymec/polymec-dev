// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
#include "geometry/sphere.h"
#include "geometry/union.h"
#include "generate_octave_script_for_surface.h"

void test_construct(void** state)
{
  // Create two spheres with inward normals.
  point_t x1 = {-0.5, 0.0, 0.0}, x2 = {0.5, 0.0, 0.0};
  sp_func_t* s1 = sphere_new(&x1, 0.25, INWARD_NORMAL);
  sp_func_t* s2 = sphere_new(&x2, 0.25, INWARD_NORMAL);

  // Now construct their union.
  sp_func_t* u = union_new2(s1, s2);
  assert_true(sp_func_num_comp(u) == 1);
  assert_true(sp_func_has_deriv(u, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  point_t x1 = {-0.5, 0.0, 0.0}, x2 = {0.5, 0.0, 0.0};
  sp_func_t* s1 = sphere_new(&x1, 0.25, INWARD_NORMAL);
  sp_func_t* s2 = sphere_new(&x2, 0.25, INWARD_NORMAL);
  sp_func_t* u = union_new2(s1, s2);
  bbox_t bbox = {.x1 = -2.0, .x2 = 2.0, .y1 = -2.0, .y2 = 2.0, .z1 = -2.0, .z2 = 2.0};
  generate_octave_script_for_surface(u, 40, &bbox, "test_union.m");
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
