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

#define DEFINE_VECTOR_TEST(vector_name, element) \
void test_##vector_name##_ctor(void** state) \
{ \
  vector_name##_t* v = vector_name##_new(10); \
  assert_int_equal(10, v->size); \
  assert_int_equal(16, v->capacity); \
  assert_true(v->data != NULL); \
  assert_true(v->arena == NULL); \
  vector_name##_free(v); \
} \
\
void test_##vector_name##_append(void** state) \
{ \
  vector_name##_t* v = vector_name##_new(0); \
  for (int i = 0; i < 10; ++i) \
    vector_name##_append(v, element##_values[i]); \
  assert_int_equal(10, v->size); \
  assert_int_equal(16, v->capacity); \
  for (int i = 0; i < 10; ++i) \
  { \
    assert_true(element##_values[i] == v->data[i]); \
  } \
  vector_name##_free(v); \
} \
\
static void count_##element(element e) \
{ \
  element##_counter++; \
} \
\
void test_##vector_name##_foreach(void** state) \
{ \
  vector_name##_t* v = vector_name##_new(0); \
  for (int i = 0; i < 10; ++i) \
    vector_name##_append(v, element##_values[i]); \
  assert_int_equal(10, v->size); \
  assert_int_equal(16, v->capacity); \
  assert_int_equal(0, element##_counter); \
  vector_name##_foreach(v, count_##element); \
  assert_int_equal(10, element##_counter); \
  vector_name##_free(v); \
} \
\
void test_##element##_vector(void** state) \
{ \
  test_##vector_name##_ctor(state); \
  test_##vector_name##_append(state); \
  test_##vector_name##_foreach(state); \
}

DEFINE_VECTOR_TEST(int_vector, int)
DEFINE_VECTOR_TEST(double_vector, double)

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
