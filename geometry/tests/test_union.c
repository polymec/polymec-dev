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
