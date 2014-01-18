// Copyright (c) 2012-2014, Jeffrey N. Johnson
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
#include "core/tuple.h"

int int_values[] = {0, 1, 2, 3, 4, 5};
int int_greater_values[] = {1, 2, 3, 4, 5, 6};
real_t real_values[] = {0., 1., 2., 3., 4., 5.};
real_t real_greater_values[] = {1., 2., 3., 4., 5., 6.};

#define DEFINE_TUPLE_TEST(tuple_name, element) \
void test_##tuple_name##_ctor(void** state) \
{ \
  tuple_name##_value_t* t = tuple_name##_new(6); \
  assert_int_equal(6, tuple_name##_length(t)); \
  for (int i = 0; i < 6; ++i) \
    t[i] = element##_values[i]; \
  assert_int_equal(6, tuple_name##_length(t)); \
  tuple_name##_free(t); \
} \
\
void test_##tuple_name##_cmp(void** state) \
{ \
  tuple_name##_value_t* t1 = tuple_name##_new(6); \
  tuple_name##_value_t* t2 = tuple_name##_new(6); \
  for (int i = 0; i < 6; ++i) \
  { \
    t1[i] = element##_values[i]; \
    t2[i] = element##_values[i]; \
  } \
  assert_int_equal(6, tuple_name##_length(t1)); \
  assert_int_equal(6, tuple_name##_length(t2)); \
  assert_int_equal(0, tuple_name##_cmp(t1, t2)); \
  for (int i = 0; i < 6; ++i) \
  { \
    t2[i] = element##_greater_values[i]; \
    assert_int_equal(-1, tuple_name##_cmp(t1, t2)); \
    t2[i] = element##_values[i]; \
    t1[i] = element##_greater_values[i]; \
    assert_int_equal(1, tuple_name##_cmp(t1, t2)); \
    t1[i] = element##_values[i]; \
  } \
  tuple_name##_free(t1); \
  tuple_name##_free(t2); \
} \
\
void test_##tuple_name##_hash(void** state) \
{ \
  tuple_name##_value_t* t1 = tuple_name##_new(6); \
  tuple_name##_value_t* t2 = tuple_name##_new(6); \
  for (int i = 0; i < 6; ++i) \
  { \
    t1[i] = element##_values[i]; \
    t2[i] = element##_greater_values[i]; \
  } \
  int h1 = tuple_name##_hash(t1); \
  int h2 = tuple_name##_hash(t2); \
  assert_int_not_equal(h1, h2); \
  tuple_name##_free(t1); \
  tuple_name##_free(t2); \
} \
void test_##tuple_name(void** state) \
{ \
  test_##tuple_name##_ctor(state); \
  test_##tuple_name##_cmp(state); \
  test_##tuple_name##_hash(state); \
}

DEFINE_TUPLE_TEST(int_tuple, int)
DEFINE_TUPLE_TEST(real_tuple, real)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_tuple),
    unit_test(test_real_tuple)
  };
  return run_tests(tests);
}
