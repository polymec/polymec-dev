// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "model/interpreter.h"

void test_interpreter_with_long_string(void** state)
{
  static const char* test_string = 
    "f = constant_function(5)\n"
    "F = constant_function(1, 2, 3)\n"
    "V = vector_function(constant_function(1), constant_function(2), constant_function(3))\n"
    "function Z(x, y, z, t)\n"
    "  return math.sqrt(x*x + y*y + z*z - t*t)\n"
    "end\n"
    "g = 2\n"
    "h = 3.0\n"
    "i = 'string cheese'\n";

  interpreter_validation_t valid_inputs[] = {{"f", INTERPRETER_SCALAR_FUNCTION, REQUIRED},
                                             {"F", INTERPRETER_VECTOR_FUNCTION, REQUIRED},
                                             {"V", INTERPRETER_VECTOR_FUNCTION, REQUIRED},
                                             {"Z", INTERPRETER_SCALAR_FUNCTION, REQUIRED},
                                             {"g", INTERPRETER_NUMBER, REQUIRED},
                                             {"h", INTERPRETER_NUMBER, REQUIRED},
                                             {"i", INTERPRETER_STRING, REQUIRED},
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
  double t = 0.0, five;
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

  assert_true(interpreter_contains(interp, "Z", INTERPRETER_SCALAR_FUNCTION));
  assert_true(interpreter_get_scalar_function(interp, "Z") != NULL);
  assert_true(!interpreter_contains(interp, "Z", INTERPRETER_VECTOR_FUNCTION));
  assert_true(!interpreter_contains(interp, "Z", INTERPRETER_NUMBER));
  assert_true(interpreter_get_number(interp, "Z") == -FLT_MAX);
  assert_true(!interpreter_contains(interp, "Z", INTERPRETER_STRING));
  assert_true(!interpreter_contains(interp, "Z", INTERPRETER_MESH));
  st_func_t* Z = interpreter_get_scalar_function(interp, "Z");
  assert_int_equal(1, st_func_num_comp(Z));
  assert_true(!st_func_is_homogeneous(Z));
  assert_true(!st_func_is_constant(Z));
  point_set(&x, 1.0, 2.0, 3.0);
  t = 2.0;
  real_t answer;
  st_func_eval(Z, &x, t, &answer);
  printf("Z = %g\n", answer);
  assert_true(fabs(answer - sqrt(10.0)) < 1e-15);

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
  interpreter_validation_t valid_inputs[] = {{"p", INTERPRETER_POINT, REQUIRED},
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
  interpreter_validation_t valid_inputs[] = {{"v", INTERPRETER_VECTOR, REQUIRED},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "v", INTERPRETER_VECTOR));
  assert_true(!interpreter_contains(interp, "w", INTERPRETER_VECTOR));
  assert_true(interpreter_get_vector(interp, "v") != NULL);
  assert_true(interpreter_get_vector(interp, "w") == NULL);
  vector_t* v = interpreter_get_vector(interp, "v");
  assert_true(fabs(v->x - 1.0) < 1e-15);
  assert_true(fabs(v->y - 2.0) < 1e-15);
  assert_true(fabs(v->z - 3.0) < 1e-15);
  interpreter_free(interp);
}

void test_boundingbox_parsing(void** state)
{
  static const char* test_string = "b = bounding_box{x1 = 0.0, x2 = 1.0, y1 = 0.0, y2 = 1.0, z1 = 0.0, z2 = 1.0}";
  interpreter_validation_t valid_inputs[] = {{"b", INTERPRETER_BOUNDING_BOX, REQUIRED},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "b", INTERPRETER_BOUNDING_BOX));
  assert_true(!interpreter_contains(interp, "c", INTERPRETER_BOUNDING_BOX));
  assert_true(interpreter_get_bbox(interp, "b") != NULL);
  assert_true(interpreter_get_bbox(interp, "c") == NULL);
  bbox_t* b = interpreter_get_bbox(interp, "b");
  assert_true(fabs(b->x1 - 0.0) < 1e-15);
  assert_true(fabs(b->x2 - 1.0) < 1e-15);
  assert_true(fabs(b->y1 - 0.0) < 1e-15);
  assert_true(fabs(b->y2 - 1.0) < 1e-15);
  assert_true(fabs(b->z1 - 0.0) < 1e-15);
  assert_true(fabs(b->z2 - 1.0) < 1e-15);
  interpreter_free(interp);
}

void test_pointlist_parsing(void** state)
{
  static const char* test_string = 
    "pts = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, "
           "{0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}";
  interpreter_validation_t valid_inputs[] = {{"pts", INTERPRETER_POINT_LIST, REQUIRED},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "pts", INTERPRETER_POINT_LIST));
  assert_true(!interpreter_contains(interp, "c", INTERPRETER_POINT_LIST));
  int num_pts;
  assert_true(interpreter_get_pointlist(interp, "pts", &num_pts) != NULL);
  assert_true(interpreter_get_pointlist(interp, "c", &num_pts) == NULL);
  point_t* pts = interpreter_get_pointlist(interp, "pts", &num_pts);
  assert_int_equal(8, num_pts);
  static point_t real_pts[8] = {{.x = 0.0, .y = 0.0, .z = 0.0},
                                {.x = 1.0, .y = 0.0, .z = 0.0},
                                {.x = 0.0, .y = 1.0, .z = 0.0},
                                {.x = 1.0, .y = 1.0, .z = 0.0},
                                {.x = 0.0, .y = 0.0, .z = 1.0},
                                {.x = 1.0, .y = 0.0, .z = 1.0},
                                {.x = 0.0, .y = 1.0, .z = 1.0},
                                {.x = 1.0, .y = 1.0, .z = 1.0}};
  for (int i = 0; i < 8; ++i)
  {
    assert_true(point_distance(&pts[i], &real_pts[i]) < 1e-15);
  }
  polymec_free(pts);
  interpreter_free(interp);
}

void test_vectorlist_parsing(void** state)
{
  static const char* test_string = 
    "vecs = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, "
            "{0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}}";
  interpreter_validation_t valid_inputs[] = {{"vecs", INTERPRETER_VECTOR_LIST, REQUIRED},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "vecs", INTERPRETER_VECTOR_LIST));
  assert_true(!interpreter_contains(interp, "c", INTERPRETER_VECTOR_LIST));
  int num_vecs;
  assert_true(interpreter_get_vectorlist(interp, "vecs", &num_vecs) != NULL);
  assert_true(interpreter_get_vectorlist(interp, "c", &num_vecs) == NULL);
  vector_t* vecs = interpreter_get_vectorlist(interp, "vecs", &num_vecs);
  assert_int_equal(8, num_vecs);
  static point_t real_vecs[8] = {{.x = 0.0, .y = 0.0, .z = 0.0},
                                 {.x = 1.0, .y = 0.0, .z = 0.0},
                                 {.x = 0.0, .y = 1.0, .z = 0.0},
                                 {.x = 1.0, .y = 1.0, .z = 0.0},
                                 {.x = 0.0, .y = 0.0, .z = 1.0},
                                 {.x = 1.0, .y = 0.0, .z = 1.0},
                                 {.x = 0.0, .y = 1.0, .z = 1.0},
                                 {.x = 1.0, .y = 1.0, .z = 1.0}};
  for (int i = 0; i < 8; ++i)
  {
    assert_true(fabs(vecs[i].x - real_vecs[i].x) < 1e-15);
    assert_true(fabs(vecs[i].y - real_vecs[i].y) < 1e-15);
    assert_true(fabs(vecs[i].z - real_vecs[i].z) < 1e-15);
  }
  polymec_free(vecs);
  interpreter_free(interp);
}

void test_table_parsing(void** state)
{
  static const char* test_string = "tab = {a = 1, b = 'bob', c = 2.0}";
  interpreter_validation_t valid_inputs[] = {{"tab", INTERPRETER_TABLE, REQUIRED},
                                             END_OF_VALID_INPUTS};
  interpreter_t* interp = interpreter_new(valid_inputs);
  interpreter_parse_string(interp, (char*)test_string);

  assert_true(interpreter_contains(interp, "tab", INTERPRETER_TABLE));
  assert_true(!interpreter_contains(interp, "c", INTERPRETER_TABLE));
  assert_true(interpreter_get_table(interp, "tab") != NULL);
  assert_true(interpreter_get_table(interp, "c") == NULL);
  string_ptr_unordered_map_t* tab = interpreter_get_table(interp, "tab");
  assert_true(tab != NULL);
  assert_int_equal(3, tab->size);
  double** a = (double**)string_ptr_unordered_map_get(tab, "a");
  assert_true(a != NULL);
  assert_true(**a == 1.0);
  char** b = (char**)string_ptr_unordered_map_get(tab, "b");
  assert_true(b != NULL);
  assert_int_equal(0, strcmp(*b, "bob"));
  double** c = (double**)string_ptr_unordered_map_get(tab, "c");
  assert_true(c != NULL);
  assert_true(**c == 2.0);
  string_ptr_unordered_map_free(tab);
  interpreter_free(interp);
}

void test_scalarfunction_parsing(void** state)
{
  point_t x = {.x = 1.0, .y = 2.0, .z = 3.0};
  real_t val, t = 1.0;

  // Constant function.
  {
    static const char* test_string = "f = 1";
    interpreter_validation_t valid_inputs[] = {{"f", INTERPRETER_SCALAR_FUNCTION, REQUIRED},
                                               END_OF_VALID_INPUTS};
    interpreter_t* interp = interpreter_new(valid_inputs);
    interpreter_parse_string(interp, (char*)test_string);

    assert_true(interpreter_contains(interp, "f", INTERPRETER_SCALAR_FUNCTION));
    st_func_t* f = interpreter_get_scalar_function(interp, "f");
    assert_true(f != NULL);
    assert_true(st_func_num_comp(f) == 1);
    st_func_eval(f, &x, t, &val);
    assert_true(fabs(val - 1.0) < 1e-15);
    interpreter_free(interp);
  }

  // Constant function (in C).
  {
    static const char* test_string = "f = constant_function(1)";
    interpreter_validation_t valid_inputs[] = {{"f", INTERPRETER_SCALAR_FUNCTION, REQUIRED},
                                               END_OF_VALID_INPUTS};
    interpreter_t* interp = interpreter_new(valid_inputs);
    interpreter_parse_string(interp, (char*)test_string);

    assert_true(interpreter_contains(interp, "f", INTERPRETER_SCALAR_FUNCTION));
    st_func_t* f = interpreter_get_scalar_function(interp, "f");
    assert_true(f != NULL);
    assert_true(st_func_num_comp(f) == 1);
    st_func_eval(f, &x, t, &val);
    assert_true(fabs(val - 1.0) < 1e-15);
    interpreter_free(interp);
  }

  // Lua function.
  {
    static const char* test_string = \
      "function f(x, y, z, t)\n"
      "  return math.sqrt(x*x + y*y + z*z - t*t)\n"
      "end\n";
    interpreter_validation_t valid_inputs[] = {{"f", INTERPRETER_SCALAR_FUNCTION, REQUIRED},
                                               END_OF_VALID_INPUTS};
    interpreter_t* interp = interpreter_new(valid_inputs);
    interpreter_parse_string(interp, (char*)test_string);

    assert_true(interpreter_contains(interp, "f", INTERPRETER_SCALAR_FUNCTION));
    st_func_t* f = interpreter_get_scalar_function(interp, "f");
    assert_true(f != NULL);
    assert_true(st_func_num_comp(f) == 1);
    st_func_eval(f, &x, t, &val);
    assert_true(fabs(val - sqrt(13.0)) < 1e-15);
    interpreter_free(interp);
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_interpreter_with_long_string),
    cmocka_unit_test(test_point_parsing),
    cmocka_unit_test(test_vector_parsing),
    cmocka_unit_test(test_boundingbox_parsing),
    cmocka_unit_test(test_pointlist_parsing),
    cmocka_unit_test(test_vectorlist_parsing),
    cmocka_unit_test(test_table_parsing),
    cmocka_unit_test(test_scalarfunction_parsing)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
