#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/table.h"

int rows[] = {0, 1, 2, 4, 8, 16};
int cols[] = {0, 10, 20, 40, 80, 160};
int int_values[] = {0, 1, 2, 3, 4, 5};
double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_TABLE_TEST(table_name, element) \
void test_##table_name##_ctor(void** state) \
{ \
  table_name##_t* t = table_name##_new(); \
  assert_int_equal(0, t->num_rows); \
  table_name##_free(t); \
} \
\
void test_##table_name##_insert(void** state) \
{ \
  table_name##_t* t = table_name##_new(); \
  for (int i = 0; i < 5; ++i) \
    table_name##_insert(t, rows[i], cols[i], element##_values[i]); \
  assert_int_equal(5, t->num_rows); \
  for (int i = 0; i < 5; ++i) \
    table_name##_insert(t, rows[i], cols[i], element##_values[i]); \
  assert_int_equal(5, t->num_rows); \
  for (int i = 0; i < 5; ++i) \
  { \
    assert_true(table_name##_contains_row(t, rows[i])); \
    assert_true(table_name##_contains(t, rows[i], cols[i])); \
    assert_true(*table_name##_get(t, rows[i], cols[i]) == element##_values[i]); \
  } \
  table_name##_clear(t); \
  assert_int_equal(0, t->num_rows); \
  table_name##_free(t); \
} \
\
void test_##element##_table(void** state) \
{ \
  test_##table_name##_ctor(state); \
  test_##table_name##_insert(state); \
}

DEFINE_TABLE_TEST(int_table, int)
DEFINE_TABLE_TEST(double_table, double)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_table),
    unit_test(test_double_table)
  };
  return run_tests(tests);
}
