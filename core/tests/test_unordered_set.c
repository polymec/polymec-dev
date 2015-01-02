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
#include "core/unordered_set.h"

int int_values[] = {0, 1, 2, 3, 4, 5};
double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_UNORDERED_SET_TEST(set_name, element) \
void test_##set_name##_ctor(void** state) \
{ \
  set_name##_t* s = set_name##_new(); \
  assert_int_equal(0, s->size); \
  set_name##_free(s); \
} \
\
void test_##set_name##_insert(void** state) \
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
void test_##element##_unordered_set(void** state) \
{ \
  test_##set_name##_ctor(state); \
  test_##set_name##_insert(state); \
}

DEFINE_UNORDERED_SET_TEST(int_unordered_set, int)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_unordered_set)
  };
  return run_tests(tests);
}
