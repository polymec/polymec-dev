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
#include "cmockery.h"
#include "geometry/sphere.h"
#include "geometry/cylinder.h"
#include "geometry/difference.h"
#include "generate_octave_script_for_surface.h"

void test_construct(void** state)
{
  // Create a sphere and a smaller cylinder.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* s = sphere_new(&origin, 0.5, INWARD_NORMAL);
  sp_func_t* c = cylinder_new(&origin, 0.25, INWARD_NORMAL);

  // With the cylinder, bore a hole through the sphere.
  sp_func_t* diff = difference_new(s, c);
  assert_true(sp_func_num_comp(diff) == 1);
  assert_true(sp_func_has_deriv(diff, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* s = sphere_new(&origin, 0.5, INWARD_NORMAL);
  sp_func_t* c = cylinder_new(&origin, 0.25, INWARD_NORMAL);
  sp_func_t* diff = difference_new(s, c);
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  generate_octave_script_for_surface(diff, 20, &bbox, "test_difference.m");
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
