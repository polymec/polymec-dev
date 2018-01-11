// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
#include "core/unordered_set.h"

static int int_values[] = {0, 1, 2, 3, 4, 5};
//static double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_UNORDERED_SET_TEST(set_name, element) \
static void test_##set_name##_ctor(void** state) \
{ \
  set_name##_t* s = set_name##_new(); \
  assert_int_equal(0, s->size); \
  set_name##_free(s); \
} \
\
static void test_##set_name##_insert(void** state) \
{ \
  set_name##_t* s = set_name##_new(); \
  for (int i = 0; i < 5; ++i) \
    set_name##_insert(s, element##_values[i]); \
  assert_int_equal(5, s->size); \
  for (int i = 0; i < 5; ++i) \
  { \
    assert_true(set_name##_contains(s, element##_values[i])); \
  } \
  set_name##_clear(s); \
  assert_int_equal(0, s->size); \
  set_name##_free(s); \
} \
\
static void test_##element##_unordered_set(void** state) \
{ \
  test_##set_name##_ctor(state); \
  test_##set_name##_insert(state); \
}

DEFINE_UNORDERED_SET_TEST(int_unordered_set, int)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_int_unordered_set)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
