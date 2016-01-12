// Copyright (c) 2012-2016, Jeffrey N. Johnson
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

void test_int_array2(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(int)*10*10);
  DECLARE_2D_ARRAY(int, a, storage, 10, 10);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      a[i][j] = 10*i + j;
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {  
      assert_int_equal(10*i+j, a[i][j]);
    }
  }
  polymec_free(storage);
} 

void test_int_array3(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(int)*10*10*10);
  DECLARE_3D_ARRAY(int, a, storage, 10, 10, 10);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        a[i][j][k] = 100*i + 10*j + k;
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {  
      for (int k = 0; k < 10; ++k)
      {  
        assert_int_equal(100*i + 10*j + k, a[i][j][k]);
      }
    }
  }
  polymec_free(storage);
} 

void test_int_array4(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(int)*10*10*10*10);
  DECLARE_4D_ARRAY(int, a, storage, 10, 10, 10, 10);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        for (int l = 0; l < 10; ++l)
          a[i][j][k][l] = 1000*i + 100*j + 10*k + l;
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {  
      for (int k = 0; k < 10; ++k)
      {  
        for (int l = 0; l < 10; ++l)
        {  
          assert_int_equal(1000*i + 100*j + 10*k + l, a[i][j][k][l]);
        }
      }
    }
  }
  polymec_free(storage);
} 

void test_int_array5(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(int)*2*2*2*2*2);
  DECLARE_5D_ARRAY(int, a, storage, 2, 2, 2, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            a[i][j][k][l][m] = 16*i + 8*j + 4*k + 2*l + m;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {  
      for (int k = 0; k < 2; ++k)
      {  
        for (int l = 0; l < 2; ++l)
        {  
          for (int m = 0; m < 2; ++m)
          {
            assert_int_equal(16*i + 8*j + 4*k + 2*l + m, a[i][j][k][l][m]);
          }
        }
      }
    }
  }
  polymec_free(storage);
} 

void test_int_array6(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(int)*2*2*2*2*2*2);
  DECLARE_6D_ARRAY(int, a, storage, 2, 2, 2, 2, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n)
              a[i][j][k][l][m][n] = 32*i + 16*j + 8*k + 4*l + 2*m + n;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {  
      for (int k = 0; k < 2; ++k)
      {  
        for (int l = 0; l < 2; ++l)
        {  
          for (int m = 0; m < 2; ++m)
          {
            for (int n = 0; n < 2; ++n)
            {
              assert_int_equal(32*i + 16*j + 8*k + 4*l + 2*m + n, a[i][j][k][l][m][n]);
            }
          }
        }
      }
    }
  }
  polymec_free(storage);
} 

#define assert_real_equal(x, y) \
  assert_true(fabs((real_t)x - y) < 1e-15)

void test_int_array7(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(int)*2*2*2*2*2*2*2);
  DECLARE_7D_ARRAY(int, a, storage, 2, 2, 2, 2, 2, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n)
              for (int p = 0; p < 2; ++p)
                a[i][j][k][l][m][n][p] = 64*i + 32*j + 16*k + 8*l + 4*m + 2*n + p;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {  
      for (int k = 0; k < 2; ++k)
      {  
        for (int l = 0; l < 2; ++l)
        {  
          for (int m = 0; m < 2; ++m)
          {
            for (int n = 0; n < 2; ++n)
            {
              for (int p = 0; p < 2; ++p)
              {
                assert_int_equal(64*i + 32*j + 16*k + 8*l + 4*m + 2*n + p, a[i][j][k][l][m][n][p]);
              }
            }
          }
        }
      }
    }
  }
  polymec_free(storage);
} 

void test_real_array2(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(real_t)*10*10);
  DECLARE_2D_ARRAY(real_t, a, storage, 10, 10);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      a[i][j] = 10*i + j;
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {  
      assert_real_equal(10*i+j, a[i][j]);
    }
  }
  polymec_free(storage);
} 

void test_real_array3(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(real_t)*10*10*10);
  DECLARE_3D_ARRAY(real_t, a, storage, 10, 10, 10);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        a[i][j][k] = 100*i + 10*j + k;
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {  
      for (int k = 0; k < 10; ++k)
      {  
        assert_real_equal(100*i + 10*j + k, a[i][j][k]);
      }
    }
  }
  polymec_free(storage);
} 

void test_real_array4(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(real_t)*10*10*10*10);
  DECLARE_4D_ARRAY(real_t, a, storage, 10, 10, 10, 10);
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      for (int k = 0; k < 10; ++k)
        for (int l = 0; l < 10; ++l)
          a[i][j][k][l] = 1000*i + 100*j + 10*k + l;
  for (int i = 0; i < 10; ++i)
  {
    for (int j = 0; j < 10; ++j)
    {  
      for (int k = 0; k < 10; ++k)
      {  
        for (int l = 0; l < 10; ++l)
        {  
          assert_real_equal(1000*i + 100*j + 10*k + l, a[i][j][k][l]);
        }
      }
    }
  }
  polymec_free(storage);
} 

void test_real_array5(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(real_t)*2*2*2*2*2);
  DECLARE_5D_ARRAY(real_t, a, storage, 2, 2, 2, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            a[i][j][k][l][m] = 16*i + 8*j + 4*k + 2*l + m;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {  
      for (int k = 0; k < 2; ++k)
      {  
        for (int l = 0; l < 2; ++l)
        {  
          for (int m = 0; m < 2; ++m)
          {
            assert_real_equal(16*i + 8*j + 4*k + 2*l + m, a[i][j][k][l][m]);
          }
        }
      }
    }
  }
  polymec_free(storage);
} 

void test_real_array6(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(real_t)*2*2*2*2*2*2);
  DECLARE_6D_ARRAY(real_t, a, storage, 2, 2, 2, 2, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n)
              a[i][j][k][l][m][n] = 32*i + 16*j + 8*k + 4*l + 2*m + n;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {  
      for (int k = 0; k < 2; ++k)
      {  
        for (int l = 0; l < 2; ++l)
        {  
          for (int m = 0; m < 2; ++m)
          {
            for (int n = 0; n < 2; ++n)
            {
              assert_real_equal(32*i + 16*j + 8*k + 4*l + 2*m + n, a[i][j][k][l][m][n]);
            }
          }
        }
      }
    }
  }
  polymec_free(storage);
} 

void test_real_array7(void** state) 
{ 
  void* storage = polymec_malloc(sizeof(real_t)*2*2*2*2*2*2*2);
  DECLARE_7D_ARRAY(real_t, a, storage, 2, 2, 2, 2, 2, 2, 2);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        for (int l = 0; l < 2; ++l)
          for (int m = 0; m < 2; ++m)
            for (int n = 0; n < 2; ++n)
              for (int p = 0; p < 2; ++p)
                a[i][j][k][l][m][n][p] = 64*i + 32*j + 16*k + 8*l + 4*m + 2*n + p;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {  
      for (int k = 0; k < 2; ++k)
      {  
        for (int l = 0; l < 2; ++l)
        {  
          for (int m = 0; m < 2; ++m)
          {
            for (int n = 0; n < 2; ++n)
            {
              for (int p = 0; p < 2; ++p)
              {
                assert_real_equal(64*i + 32*j + 16*k + 8*l + 4*m + 2*n + p, a[i][j][k][l][m][n][p]);
              }
            }
          }
        }
      }
    }
  }
  polymec_free(storage);
} 

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
    cmocka_unit_test(test_real_array2),
    cmocka_unit_test(test_real_array3),
    cmocka_unit_test(test_real_array4),
    cmocka_unit_test(test_real_array5),
    cmocka_unit_test(test_real_array6),
    cmocka_unit_test(test_real_array7)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
