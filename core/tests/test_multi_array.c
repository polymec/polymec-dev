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
