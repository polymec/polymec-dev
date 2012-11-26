#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "interpreter.h"

const char* test_string = 
"f = constant_function(5)\n"
"g = 2\n"
"h = 3.0\n"
"i = 'string cheese'\n";

void test_interpreter(void** state)
{
  interpreter_validation_t valid_inputs[] = {{"f", INTERPRETER_FUNCTION},
                                             {"g", INTERPRETER_NUMBER},
                                             {"h", INTERPRETER_NUMBER},
                                             {"i", INTERPRETER_STRING},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "f", INTERPRETER_FUNCTION));
  assert_true(interpreter_get_function(interp, "f") != NULL);
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_NUMBER));
  assert_true(interpreter_get_number(interp, "f") == -FLT_MAX);
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_STRING));
  assert_true(!interpreter_contains(interp, "f", INTERPRETER_MESH));
  st_func_t* F = interpreter_get_function(interp, "f");
  assert_true(st_func_is_homogeneous(F));
  assert_true(st_func_is_constant(F));
  point_t x;
  double t, five;
  st_func_eval(F, &x, t, &five);
  assert_true(fabs(five - 5.0) < 1e-15);

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
  arbi_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_interpreter)
  };
  return run_tests(tests);
}
