#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/vector.h"

static const int int_values[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
static const double double_values[] = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
static int int_counter = 0;
static int double_counter = 0;

#define DEFINE_VECTOR_TEST(element) \
void test_##element##_vector_ctor(void** state) \
{ \
  element##_vector_t* v = element##_vector_new(10); \
  assert_int_equal(10, v->size); \
  assert_int_equal(16, v->capacity); \
  assert_true(v->data != NULL); \
  assert_true(v->arena != NULL); \
  element##_vector_free(v); \
} \
\
void test_##element##_vector_append(void** state) \
{ \
  element##_vector_t* v = element##_vector_new(0); \
  for (int i = 0; i < 10; ++i) \
    element##_vector_append(v, element##_values[i]); \
  assert_int_equal(10, v->size); \
  assert_int_equal(16, v->capacity); \
  for (int i = 0; i < 10; ++i) \
  { \
    assert_true(element##_values[i] == v->data[i]); \
  } \
  element##_vector_free(v); \
} \
\
static void count_##element(element e) \
{ \
  element##_counter++; \
} \
\
void test_##element##_vector_foreach(void** state) \
{ \
  element##_vector_t* v = element##_vector_new(0); \
  for (int i = 0; i < 10; ++i) \
    element##_vector_append(v, element##_values[i]); \
  assert_int_equal(10, v->size); \
  assert_int_equal(16, v->capacity); \
  assert_int_equal(0, element##_counter); \
  element##_vector_foreach(v, count_##element); \
  assert_int_equal(10, element##_counter); \
  element##_vector_free(v); \
} \
\
void test_##element##_vector(void** state) \
{ \
  test_##element##_vector_ctor(state); \
  test_##element##_vector_append(state); \
  test_##element##_vector_foreach(state); \
}

DEFINE_VECTOR_TEST(int)
DEFINE_VECTOR_TEST(double)

int main(int argc, char* argv[]) 
{
  arbi_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_vector),
    unit_test(test_double_vector)
  };
  return run_tests(tests);
}
