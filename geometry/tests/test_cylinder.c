// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "geometry/cylinder.h"
#include "generate_octave_script_for_surface.h"

void test_construct(void** state)
{
  // Create spheres with inward/outward normals.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* c1 = cylinder_new(&origin, 1.0, INWARD_NORMAL);
  assert_true(sp_func_num_comp(c1) == 1);
  assert_true(sp_func_has_deriv(c1, 1));
  sp_func_t* c2 = cylinder_new(&origin, 1.0, OUTWARD_NORMAL);
  assert_true(sp_func_num_comp(c2) == 1);
  assert_true(sp_func_has_deriv(c2, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* c = cylinder_new(&origin, 0.5, OUTWARD_NORMAL);
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  generate_octave_script_for_surface(c, 20, &bbox, "test_cylinder.m");
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
