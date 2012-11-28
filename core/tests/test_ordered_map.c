#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/ordered_map.h"

const char* keys[] = {"0", "1", "2", "3", "4", "5"};
int int_values[] = {0, 1, 2, 3, 4, 5};
double double_values[] = {0., 1., 2., 3., 4., 5.};

// Key-value destructors.
#define DEFINE_ORDERED_MAP_TEST(map_name, element) \
DEFINE_ORDERED_MAP(map_name, char*, element, strcmp) \
void test_##map_name##_ctor(void** state) \
{ \
  map_name##_t* s = map_name##_new(); \
  assert_int_equal(0, s->size); \
  map_name##_free(s); \
} \
\
void test_##map_name##_insert(void** state) \
{ \
  map_name##_t* s = map_name##_new(); \
  for (int i = 0; i < 5; ++i) \
    map_name##_insert(s, (char*)keys[i], element##_values[i]); \
  assert_int_equal(5, s->size); \
  for (int i = 0; i < 5; ++i) \
  { \
    assert_true(map_name##_find(s, (char*)keys[i]) != NULL); \
    assert_true(map_name##_value(s, (char*)keys[i]) == element##_values[i]); \
  } \
  map_name##_clear(s); \
  assert_int_equal(0, s->size); \
  map_name##_free(s); \
} \
\
void test_##element##_ordered_map(void** state) \
{ \
  test_##map_name##_ctor(state); \
  test_##map_name##_insert(state); \
}

DEFINE_ORDERED_MAP_TEST(int_ordered_map, int)
DEFINE_ORDERED_MAP_TEST(double_ordered_map, double)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_ordered_map),
    unit_test(test_double_ordered_map)
  };
  return run_tests(tests);
}
