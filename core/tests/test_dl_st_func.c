#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "cmockery.h"
#include "dl_st_func.h"

static const char* source = 
"void eval(void* context, point_t* x, double t, double* result)\n"
"{\n"
"  result[0] = (x->x + x->y + x->z) * sin(t);\n"
"  result[1] = (x->x + x->y + x->z) * cos(t);\n"
"}\n"
"bool homogeneous = false;\n"
"int num_comp = 2;\n";

void test_construction(void** state)
{
  UNUSED_ARG(state);
  
  // Compile an object file from the source.
  dl_st_func_set_compiler("cc", "-shared");
  dl_st_func_register("example", source);

  // Create an st_func and make sure it has the right metadata.
  st_func_t* func = dl_st_func_new("example");
  assert_true(func != NULL);
  assert_false(st_func_is_homogeneous(func));
  assert_false(st_func_is_constant(func));
  assert_int_equal(2, st_func_num_comp(func));

  // Try calling it.
  point_t x = {.x = 1.0, .y = 2.0, .z = 3.0};
  double result[2];
  st_func_eval(func, &x, 1.0, result);
  assert_true(fabs(result[0] - 6.0*sin(1.0)) < 1e-14);
  assert_true(fabs(result[1] - 6.0*cos(1.0)) < 1e-14);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_construction)
  };
  return run_tests(tests);
}
