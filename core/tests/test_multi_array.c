// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/multi_array.h"

void test_int_array2(void** state) 
{ 
  int** a = int_array2_new(10, 10); 
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      a[i][j] = 10*i + j;
  for (int i = 0; i < 100; ++i)
    assert_int_equal(i, a[0][i]);
  int_array2_free(a, 10, 10); 
} 

void test_int_array3(void** state) 
{ 
  int*** a = int_array3_new(10, 10, 10); 
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        a[i][j][k] = 100*i + 10*j + k;
  for (int i = 0; i < 1000; ++i)
    assert_int_equal(i, a[0][0][i]);
  int_array3_free(a, 10, 10, 10); 
} 

void test_int_array4(void** state) 
{ 
  int**** a = int_array4_new(10, 10, 10, 10); 
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        for (int l = 0; l < 10; ++l)
          a[i][j][k][l] = 1000*i + 100*j + 10*k + l;
  for (int i = 0; i < 10000; ++i)
    assert_int_equal(i, a[0][0][0][i]);
  int_array4_free(a, 10, 10, 10, 10); 
} 

void test_real_array2(void** state) 
{ 
  real_t** a = real_array2_new(10, 10); 
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      a[i][j] = (real_t)1.0*(10*i + j);
  for (int i = 0; i < 100; ++i)
    assert_true(a[0][i] == (real_t)(i));
  real_array2_free(a, 10, 10); 
} 

void test_real_array3(void** state) 
{ 
  real_t*** a = real_array3_new(10, 10, 10); 
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        a[i][j][k] = (real_t)1.0*(100*i + 10*j + k);
  for (int i = 0; i < 1000; ++i)
    assert_true(a[0][0][i] == (real_t)i);
  real_array3_free(a, 10, 10, 10); 
} 

void test_real_array4(void** state) 
{ 
  real_t**** a = real_array4_new(10, 10, 10, 10); 
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        for (int l = 0; l < 10; ++l)
          a[i][j][k][l] = (real_t)(1000*i + 100*j + 10*k + l);
  for (int i = 0; i < 10000; ++i)
    assert_true(a[0][0][0][i] == (real_t)i);
  real_array4_free(a, 10, 10, 10, 10); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_array2),
    unit_test(test_int_array3),
    unit_test(test_int_array4),
    unit_test(test_real_array2),
    unit_test(test_real_array3),
    unit_test(test_real_array4)
  };
  return run_tests(tests);
}
