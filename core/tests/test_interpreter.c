#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "interpreter.h"

const char* test_string = 
"f = constant_function(5)\n"
"F = constant_function(1, 2, 3)\n"
"V = vector_function(constant_function(1), constant_function(2), constant_function(3))\n"
"g = 2\n"
"h = 3.0\n"
"i = 'string cheese'\n";

void test_interpreter(void** state)
{
  interpreter_validation_t valid_inputs[] = {{"f", INTERPRETER_SCALAR_FUNCTION},
                                             {"F", INTERPRETER_VECTOR_FUNCTION},
                                             {"V", INTERPRETER_VECTOR_FUNCTION},
                                             {"g", INTERPRETER_NUMBER},
                                             {"h", INTERPRETER_NUMBER},
                                             {"i", INTERPRETER_STRING},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "f", INTERPRETER_SCALAR_FUNCTION));
  assert_true(interpreter_get_scalar_function(interp, "f") != NULL);
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_VECTOR_FUNCTION));
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_NUMBER));
  assert_true(interpreter_get_number(interp, "f") == -FLT_MAX);
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_STRING));
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_MESH));
  st_func_t* f = interpreter_get_scalar_function(interp, "f");
  assert_int_equal(1, st_func_num_comp(f));
  assert_true(st_func_is_homogeneous(f));
  assert_true(st_func_is_constant(f));
  point_t x;
  double t, five;
  st_func_eval(f, &x, t, &five);
  assert_true(fabs(five - 5.0) < 1e-15);

  assert_true(interpreter_contains(interp, "F", INTERPRETER_VECTOR_FUNCTION));
  assert_true(interpreter_get_vector_function(interp, "F") != NULL);
  assert_true(!interpreter_contains(interp, "F", INTERPRETER_SCALAR_FUNCTION));
  assert_true(!interpreter_contains(interp, "F", INTERPRETER_NUMBER));
  assert_true(interpreter_get_number(interp, "F") == -FLT_MAX);
  assert_true(!interpreter_contains(interp, "F", INTERPRETER_STRING));
  assert_true(!interpreter_contains(interp, "F", INTERPRETER_MESH));
  st_func_t* F = interpreter_get_vector_function(interp, "F");
  assert_int_equal(3, st_func_num_comp(F));
  assert_true(st_func_is_homogeneous(F));
  assert_true(st_func_is_constant(F));
  double one_two_three[3];
  st_func_eval(F, &x, t, one_two_three);
  assert_true(fabs(one_two_three[0] - 1.0) < 1e-15);
  assert_true(fabs(one_two_three[1] - 2.0) < 1e-15);
  assert_true(fabs(one_two_three[2] - 3.0) < 1e-15);

  assert_true(interpreter_contains(interp, "V", INTERPRETER_VECTOR_FUNCTION));
  assert_true(interpreter_get_vector_function(interp, "V") != NULL);
  st_func_t* V = interpreter_get_vector_function(interp, "V");
  assert_int_equal(3, st_func_num_comp(V));
  assert_true(st_func_is_homogeneous(V));
  assert_true(st_func_is_constant(V));
  st_func_eval(V, &x, t, one_two_three);
  assert_true(fabs(one_two_three[0] - 1.0) < 1e-15);
  assert_true(fabs(one_two_three[1] - 2.0) < 1e-15);
  assert_true(fabs(one_two_three[2] - 3.0) < 1e-15);

  assert_true(interpreter_contains(interp, "g", INTERPRETER_NUMBER));
  assert_true((int)interpreter_get_number(interp, "g") == 2);

  assert_true(interpreter_contains(interp, "h", INTERPRETER_NUMBER));
  assert_true(interpreter_get_number(interp, "h") == 3.0);

  assert_true(interpreter_contains(interp, "i", INTERPRETER_STRING));
  assert_true(interpreter_get_number(interp, "i") == -FLT_MAX);
  assert_true(!strcmp(interpreter_get_string(interp, "i"), "string cheese"));

  interpreter_free(interp);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_interpreter)
  };
  return run_tests(tests);
}
