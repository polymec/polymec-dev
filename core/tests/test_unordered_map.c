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
#include "core/unordered_map.h"
#include "core/hash_functions.h"
#include "core/comparators.h"

static const char* keys[] = {"0", "1", "2", "3", "4", "5"};
static int int_values[] = {0, 1, 2, 3, 4, 5};
static double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_UNORDERED_MAP_TEST(map_name, element) \
DEFINE_UNORDERED_MAP(map_name, char*, element, string_hash, string_equals) \
static void test_##map_name##_ctor(void** state) \
{ \
  map_name##_t* m = map_name##_new(); \
  assert_int_equal(0, m->size); \
  map_name##_free(m); \
} \
\
static void test_##map_name##_insert(void** state) \
{ \
  map_name##_t* m = map_name##_new(); \
  for (int i = 0; i < 5; ++i) \
    map_name##_insert(m, (char*)keys[i], element##_values[i]); \
  assert_int_equal(5, m->size); \
  for (int i = 0; i < 5; ++i) \
    map_name##_insert(m, (char*)keys[i], element##_values[i]); \
  assert_int_equal(5, m->size); \
  for (int i = 0; i < 5; ++i) \
  { \
    assert_true(map_name##_contains(m, (char*)keys[i])); \
    assert_true(*map_name##_get(m, (char*)keys[i]) == element##_values[i]); \
  } \
  map_name##_clear(m); \
  assert_int_equal(0, m->size); \
  map_name##_free(m); \
} \
\
static void test_##map_name##_swap(void** state) \
{ \
  map_name##_t* m = map_name##_new(); \
  for (int i = 0; i < 2; ++i) \
    map_name##_insert(m, (char*)keys[i], element##_values[i]); \
  map_name##_swap(m, (char*)keys[0], (char*)keys[1]); \
  assert_true(map_name##_contains(m, (char*)keys[0])); \
  assert_true(map_name##_contains(m, (char*)keys[1])); \
  assert_true(*map_name##_get(m, (char*)keys[0]) == element##_values[1]); \
  assert_true(*map_name##_get(m, (char*)keys[1]) == element##_values[0]); \
  map_name##_free(m); \
} \
\
static void test_##element##_unordered_map(void** state) \
{ \
  test_##map_name##_ctor(state); \
  test_##map_name##_insert(state); \
  test_##map_name##_swap(state); \
}

DEFINE_UNORDERED_MAP_TEST(int_unordered_map, int)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
DEFINE_UNORDERED_MAP_TEST(double_unordered_map, double)
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_int_unordered_map),
    cmocka_unit_test(test_double_unordered_map)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
