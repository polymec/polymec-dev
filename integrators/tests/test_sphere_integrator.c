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
#include "core/polymec.h"
#include "core/constant_st_func.h"
#include "integrators/sphere_integrator.h"

void test_ctor(void** state)
{
  sphere_integrator_rule_t rules[] = { GAUSS_LEGENDRE, GAUSS_RADAU, GAUSS_LOBATTO};
  for (int order = 1; order <= 4; ++order)
  {
    for (int rule = 0; rule < 3; ++rule)
    {
      sphere_integrator_t* I = sphere_integrator_new(order, rules[rule]);
      sphere_integrator_free(I);
    }
  }
}

void test_cap_area(void** state)
{
  double o = 1.0;
  sp_func_t* one = constant_sp_func_new(1, &o);
  vector_t e3 = {.x = 0.0, .y = 0.0, .z = 1.0};
  sphere_integrator_rule_t rules[] = { GAUSS_LEGENDRE, GAUSS_RADAU, GAUSS_LOBATTO};
  for (int order = 1; order <= 4; ++order)
  {
    for (int rule = 0; rule < 3; ++rule)
    {
      sphere_integrator_t* I = sphere_integrator_new(order, rules[rule]);

      static const double radii[] = {1.0, 2.0, 3.0};
      for (int r = 0; r < 3; ++r)
      {
        double radius = radii[r];

        // This should give zero.
        double area;
        sphere_integrator_cap(I, radius, one, &e3, 0.0, &area);
        assert_true(fabs(area) < 1e-8);

        // This should give the area of the entire sphere.
        sphere_integrator_cap(I, radius, one, &e3, M_PI, &area);
        printf("%g - %g = %g\n", area, 4.0*M_PI*radius*radius, area - 4.0*M_PI*radius*radius);
        assert_true(fabs(area - 4.0*M_PI*radius*radius) < 1e-8);

        // This should give half the area of the sphere.
        sphere_integrator_cap(I, radius, one, &e3, 0.5*M_PI, &area);
        assert_true(fabs(area - 2.0*M_PI*radius*radius) < 1e-8);
      }

      sphere_integrator_free(I);
    }
  }
  one = NULL;
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_cap_area)
//    unit_test(test_cap_at_time)
  };
  return run_tests(tests);
}
