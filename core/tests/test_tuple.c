#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/tuple.h"

int int_values[] = {0, 1, 2, 3, 4, 5};
int int_greater_values[] = {1, 2, 3, 4, 5, 6};
double double_values[] = {0., 1., 2., 3., 4., 5.};
double double_greater_values[] = {1., 2., 3., 4., 5., 6.};

#define DEFINE_TUPLE_TEST(tuple_name, element) \
void test_##tuple_name##_ctor(void** state) \
{ \
  element* t = tuple_name##_new(6); \
  assert_int_equal(6, tuple_name##_length(t)); \
  for (int i = 0; i < 6; ++i) \
    t[i] = element##_values[i]; \
  assert_int_equal(6, tuple_name##_length(t)); \
  tuple_name##_free(t); \
} \
\
void test_##tuple_name##_cmp(void** state) \
{ \
  element* t1 = tuple_name##_new(6); \
  element* t2 = tuple_name##_new(6); \
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
  element* t1 = tuple_name##_new(6); \
  element* t2 = tuple_name##_new(6); \
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
DEFINE_TUPLE_TEST(double_tuple, double)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_tuple),
    unit_test(test_double_tuple)
  };
  return run_tests(tests);
}
