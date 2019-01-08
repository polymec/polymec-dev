// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/table.h"

static int rows[] = {0, 1, 2, 4, 8, 16};
static int cols[] = {0, 10, 20, 40, 80, 160};
static int int_values[] = {0, 1, 2, 3, 4, 5};
static real_t real_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_TABLE_TEST(table_name, element) \
static void test_##table_name##_ctor(void** state) \
{ \
  table_name##_t* t = table_name##_new(); \
  assert_int_equal(0, t->num_rows); \
  table_name##_free(t); \
} \
\
static void test_##table_name##_insert(void** state) \
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
static void test_##table_name##_tridiag(void** state) \
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
static void test_##element##_table(void** state) \
{ \
  test_##table_name##_ctor(state); \
  test_##table_name##_insert(state); \
  test_##table_name##_tridiag(state); \
}

DEFINE_TABLE_TEST(int_table, int)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
DEFINE_TABLE_TEST(real_table, real)
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_int_table),
    cmocka_unit_test(test_real_table)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
