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
#include "core/hash_map.h"
#include "core/hash_functions.h"
#include "core/comparators.h"

const char* keys[] = {"0", "1", "2", "3", "4", "5"};
int int_values[] = {0, 1, 2, 3, 4, 5};
double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_HASH_MAP_TEST(map_name, element) \
DEFINE_HASH_MAP(map_name, char*, element, string_hash, string_equals) \
void test_##map_name##_ctor(void** state) \
{ \
  map_name##_t* m = map_name##_new(); \
  assert_int_equal(0, m->size); \
  map_name##_free(m); \
} \
\
void test_##map_name##_insert(void** state) \
{ \
  map_name##_t* m = map_name##_new(); \
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
void test_##element##_hash_map(void** state) \
{ \
  test_##map_name##_ctor(state); \
  test_##map_name##_insert(state); \
}

DEFINE_HASH_MAP_TEST(int_hash_map, int)
DEFINE_HASH_MAP_TEST(double_hash_map, double)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_hash_map),
    unit_test(test_double_hash_map)
  };
  return run_tests(tests);
}
