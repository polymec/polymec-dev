// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
#include "core/declare_nd_array.h"

#define DEFINE_ND_ARRAY_TESTS(X) \
static void test_##X##_array2(void** state) \
{ \
  void* storage = polymec_malloc(sizeof(X)*9*10); \
  for (int i = 0; i < 9*10; ++i) \
    ((X*)storage)[i] = (X)i; \
  DECLARE_2D_ARRAY(X, a, storage, 9, 10); \
  for (int i = 0; i < 9; ++i) \
  { \
    for (int j = 0; j < 10; ++j) \
    { \
      assert_##X##_equal(a[i][j], ((X*)storage)[ARRAY_INDEX_2D(9, 10, i, j)]); \
    } \
  } \
  polymec_free(storage); \
} \
\
static void test_##X##_array3(void** state) \
{ \
  void* storage = polymec_malloc(sizeof(X)*8*9*10); \
  for (int i = 0; i < 8*9*10; ++i) \
    ((X*)storage)[i] = (X)i; \
  DECLARE_3D_ARRAY(X, a, storage, 8, 9, 10); \
  for (int i = 0; i < 8; ++i) \
  { \
    for (int j = 0; j < 9; ++j) \
    { \
      for (int k = 0; k < 10; ++k) \
      { \
        assert_##X##_equal(a[i][j][k], ARRAY_INDEX_3D(8, 9, 10, i, j, k)); \
      } \
    } \
  } \
  polymec_free(storage); \
} \
\
static void test_##X##_array4(void** state) \
{ \
  void* storage = polymec_malloc(sizeof(X)*7*8*9*10); \
  for (int i = 0; i < 7*8*9*10; ++i) \
    ((X*)storage)[i] = (X)i; \
  DECLARE_4D_ARRAY(X, a, storage, 7, 8, 9, 10); \
  for (int i = 0; i < 7; ++i) \
  { \
    for (int j = 0; j < 8; ++j) \
    { \
      for (int k = 0; k < 9; ++k) \
      { \
        for (int l = 0; l < 10; ++l) \
        { \
          assert_##X##_equal(a[i][j][k][l], ARRAY_INDEX_4D(7, 8, 9, 10, i, j, k, l)); \
        } \
      } \
    } \
  } \
  polymec_free(storage); \
} \
\
static void test_##X##_array5(void** state) \
{ \
  void* storage = polymec_malloc(sizeof(X)*1*2*3*4*5); \
  for (int i = 0; i < 1*2*3*4*5; ++i) \
    ((X*)storage)[i] = i; \
  DECLARE_5D_ARRAY(X, a, storage, 1, 2, 3, 4, 5); \
  for (int i = 0; i < 1; ++i) \
  { \
    for (int j = 0; j < 2; ++j) \
    { \
      for (int k = 0; k < 3; ++k) \
      { \
        for (int l = 0; l < 4; ++l) \
        { \
          for (int m = 0; m < 5; ++m) \
          { \
            assert_##X##_equal(a[i][j][k][l][m], ARRAY_INDEX_5D(1, 2, 3, 4, 5, i, j, k, l, m)); \
          } \
        } \
      } \
    } \
  } \
  polymec_free(storage); \
} \
\
static void test_##X##_array6(void** state) \
{ \
  void* storage = polymec_malloc(sizeof(X)*1*2*3*4*5*6); \
  for (int i = 0; i < 1*2*3*4*5*6; ++i) \
    ((X*)storage)[i] = i; \
  DECLARE_6D_ARRAY(X, a, storage, 1, 2, 3, 4, 5, 6); \
  for (int i = 0; i < 1; ++i) \
  { \
    for (int j = 0; j < 2; ++j) \
    { \
      for (int k = 0; k < 3; ++k) \
      { \
        for (int l = 0; l < 4; ++l) \
        { \
          for (int m = 0; m < 5; ++m) \
          { \
            for (int n = 0; n < 6; ++n) \
            { \
              assert_##X##_equal(a[i][j][k][l][m][n], ARRAY_INDEX_6D(1, 2, 3, 4, 5, 6, i, j, k, l, m, n)); \
            } \
          } \
        } \
      } \
    } \
  } \
  polymec_free(storage); \
} \
\
static void test_##X##_array7(void** state) \
{ \
  void* storage = polymec_malloc(sizeof(X)*1*2*3*4*5*6*7); \
  for (int i = 0; i < 1*2*3*4*5*6*7; ++i) \
    ((X*)storage)[i] = i; \
  DECLARE_7D_ARRAY(X, a, storage, 1, 2, 3, 4, 5, 6, 7); \
  for (int i = 0; i < 1; ++i) \
  { \
    for (int j = 0; j < 2; ++j) \
    { \
      for (int k = 0; k < 3; ++k) \
      { \
        for (int l = 0; l < 4; ++l) \
        { \
          for (int m = 0; m < 5; ++m) \
          { \
            for (int n = 0; n < 6; ++n) \
            { \
              for (int p = 0; p < 7; ++p) \
              { \
                assert_##X##_equal(a[i][j][k][l][m][n][p], ARRAY_INDEX_7D(1, 2, 3, 4, 5, 6, 7, i, j, k, l, m, n, p)); \
              } \
            } \
          } \
        } \
      } \
    } \
  } \
  polymec_free(storage); \
} \

DEFINE_ND_ARRAY_TESTS(int)
#define assert_real_t_equal(x, y) \
  assert_true(reals_nearly_equal((real_t)x, y, 1e-15))
DEFINE_ND_ARRAY_TESTS(real_t)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_int_array2),
    cmocka_unit_test(test_int_array3),
    cmocka_unit_test(test_int_array4),
    cmocka_unit_test(test_int_array5),
    cmocka_unit_test(test_int_array6),
    cmocka_unit_test(test_int_array7),
    cmocka_unit_test(test_real_t_array2),
    cmocka_unit_test(test_real_t_array3),
    cmocka_unit_test(test_real_t_array4),
    cmocka_unit_test(test_real_t_array5),
    cmocka_unit_test(test_real_t_array6),
    cmocka_unit_test(test_real_t_array7)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
