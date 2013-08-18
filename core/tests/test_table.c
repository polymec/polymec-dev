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
void test_##table_name##_tridiag(void** state) \
{ \
  table_name##_t* t = table_name##_new(); \
  for (int i = 0; i < 32; ++i) \
  { \
    table_name##_insert(t, i, i, (table_name##_value_t)i); \
    table_name##_insert(t, i, i-1, (table_name##_value_t)(i-1)); \
    table_name##_insert(t, i, i+1, (table_name##_value_t)(i+1)); \
  } \
  assert_int_equal(32, t->num_rows); \
  { \
    table_name##_cell_pos_t pos = table_name##_start(t); \
    int i, j, nnz = 0; \
    table_name##_value_t tij; \
    while (table_name##_next_cell(t, &pos, &i, &j, &tij)) \
    { \
      assert_true((i == j) || (i == (j+1)) || (i == (j-1))); \
      assert_true((table_name##_value_t)tij == j); \
      ++nnz; \
    } \
    assert_int_equal(32*3, nnz); \
  } \
  for (int i = 0; i < 32; ++i) \
  { \
    table_name##_insert(t, i, i, (table_name##_value_t)i); \
    table_name##_insert(t, i, i-1, (table_name##_value_t)(i-1)); \
    table_name##_insert(t, i, i+1, (table_name##_value_t)(i+1)); \
  } \
  assert_int_equal(32, t->num_rows); \
  { \
    table_name##_cell_pos_t pos = table_name##_start(t); \
    int i, j, nnz = 0; \
    table_name##_value_t tij; \
    while (table_name##_next_cell(t, &pos, &i, &j, &tij)) \
    { \
      assert_true((i == j) || (i == (j+1)) || (i == (j-1))); \
      assert_true((table_name##_value_t)tij == j); \
      ++nnz; \
    } \
    assert_int_equal(32*3, nnz); \
  } \
  table_name##_clear(t); \
  assert_int_equal(0, t->num_rows); \
  { \
    table_name##_cell_pos_t pos = table_name##_start(t); \
    int i, j, nnz = 0; \
    table_name##_value_t tij; \
    while (table_name##_next_cell(t, &pos, &i, &j, &tij)) \
      ++nnz; \
    assert_int_equal(0, nnz); \
  } \
  table_name##_free(t); \
} \
\
void test_##element##_table(void** state) \
{ \
  test_##table_name##_ctor(state); \
  test_##table_name##_insert(state); \
  test_##table_name##_tridiag(state); \
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
