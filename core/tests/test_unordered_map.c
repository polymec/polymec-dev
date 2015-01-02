// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/unordered_map.h"
#include "core/hash_functions.h"
#include "core/comparators.h"

const char* keys[] = {"0", "1", "2", "3", "4", "5"};
int int_values[] = {0, 1, 2, 3, 4, 5};
double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_UNORDERED_MAP_TEST(map_name, element) \
DEFINE_UNORDERED_MAP(map_name, char*, element, string_hash, string_equals) \
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
void test_##element##_unordered_map(void** state) \
{ \
  test_##map_name##_ctor(state); \
  test_##map_name##_insert(state); \
}

DEFINE_UNORDERED_MAP_TEST(int_unordered_map, int)
DEFINE_UNORDERED_MAP_TEST(double_unordered_map, double)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_unordered_map),
    unit_test(test_double_unordered_map)
  };
  return run_tests(tests);
}
