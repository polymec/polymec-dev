#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/interpreter.h"

void test_interpreter_with_long_string(void** state)
{
  static const char* test_string = 
    "f = constant_function(5)\n"
    "F = constant_function(1, 2, 3)\n"
    "V = vector_function(constant_function(1), constant_function(2), constant_function(3))\n"
    "g = 2\n"
    "h = 3.0\n"
    "i = 'string cheese'\n";

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

void test_point_parsing(void** state)
{
  static const char* test_string = "p = {1.0, 2.0, 3.0}";
  interpreter_validation_t valid_inputs[] = {{"p", INTERPRETER_POINT},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "p", INTERPRETER_POINT));
  assert_true(!interpreter_contains(interp, "q", INTERPRETER_POINT));
  assert_true(interpreter_get_point(interp, "p") != NULL);
  assert_true(interpreter_get_point(interp, "q") == NULL);
  point_t* p = interpreter_get_point(interp, "p");
  assert_true(fabs(p->x - 1.0) < 1e-15);
  assert_true(fabs(p->y - 2.0) < 1e-15);
  assert_true(fabs(p->z - 3.0) < 1e-15);
  interpreter_free(interp);
}

void test_vector_parsing(void** state)
{
  static const char* test_string = "v = {1.0, 2.0, 3.0}";
  interpreter_validation_t valid_inputs[] = {{"v", INTERPRETER_VECTOR},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "v", INTERPRETER_VECTOR));
  assert_true(!interpreter_contains(interp, "w", INTERPRETER_POINT));
  assert_true(interpreter_get_vector(interp, "v") != NULL);
  assert_true(interpreter_get_vector(interp, "w") == NULL);
  vector_t* v = interpreter_get_vector(interp, "v");
  assert_true(fabs(v->x - 1.0) < 1e-15);
  assert_true(fabs(v->y - 2.0) < 1e-15);
  assert_true(fabs(v->z - 3.0) < 1e-15);
  interpreter_free(interp);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_interpreter_with_long_string),
    unit_test(test_point_parsing),
    unit_test(test_vector_parsing)
  };
  return run_tests(tests);
}
