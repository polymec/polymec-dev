#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "geometry/sphere.h"
#include "generate_octave_script_for_surface.h"

void test_construct(void** state)
{
  // Create spheres with inward/outward normals.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* s1 = sphere_new(&origin, 1.0, INWARD_NORMAL);
  assert_true(sp_func_num_comp(s1) == 1);
  assert_true(sp_func_has_deriv(s1, 1));
  sp_func_t* s2 = sphere_new(&origin, 1.0, OUTWARD_NORMAL);
  assert_true(sp_func_num_comp(s2) == 1);
  assert_true(sp_func_has_deriv(s2, 1));
}

void test_plot(void** state)
{
  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  point_t origin = {0.0, 0.0, 0.0};
  sp_func_t* s = sphere_new(&origin, 0.5, OUTWARD_NORMAL);
  bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};
  generate_octave_script_for_surface(s, 20, &bbox, "test_sphere.m");
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
