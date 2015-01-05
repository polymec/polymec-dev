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
#include "core/polymec.h"
#include "core/special_functions.h"

void test_find_jn_roots(void** state)
{
  static double jn_roots[6][5] = 
    {{2.4048, 5.5201, 8.6537, 11.7915, 14.9309},    // J0
     {3.8317, 7.0156, 10.1735, 13.3237, 16.4706},   // J1
     {5.1356, 8.4172, 11.6198, 14.7960, 17.9598},   // J2
     {6.3802, 9.7610, 13.0152, 16.2235, 19.4094},   // J3
     {7.5883, 11.0647, 14.3725, 17.6160, 20.8269},  // J4
     {8.7715, 12.3386, 15.7002, 18.9801, 22.2178}}; // J5
   
  for (int n = 0; n <= 5; ++n)
  {
    double roots[5];
    find_jn_roots(n, 5, roots);
    for (int k = 0; k < 5; ++k)
    {
      assert_true(fabs(roots[k] - jn_roots[n][k])/jn_roots[n][k] < 1e-4);
    }
  }
}

void test_find_yn_roots(void** state)
{
  // These reference roots were taken from scipy.special.yn_zeros.
  static double yn_roots[6][5] = 
    {{0.8936,   3.9577,   7.0861,  10.2223,  13.3611},  // Y0
     {2.1971,   5.4297,   8.5960,  11.7492,  14.8974},  // Y1
     {3.3842,   6.7938,  10.0235,  13.2100,  16.3790},  // Y2
     {4.5270,   8.0976,  11.3965,  14.6231,  17.8185},  // Y3
     {5.6451,   9.3616,  12.7301,  15.9996,  19.2244},  // Y4
     {6.7472,  10.5972,  14.0338,  17.3471,  20.6029}}; // Y5

  for (int n = 0; n <= 5; ++n)
  {
    double roots[5];
    find_yn_roots(n, 5, roots);
    for (int k = 0; k < 5; ++k)
    {
      assert_true(fabs(roots[k] - yn_roots[n][k])/yn_roots[n][k] < 1e-4);
    }
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_find_jn_roots),
//    unit_test(test_j0),
//    unit_test(test_j1),
//    unit_test(test_jn),
//    unit_test(test_dj0dx),
//    unit_test(test_dj1dx),
//    unit_test(test_djndx),
    unit_test(test_find_yn_roots),
//    unit_test(test_y0),
//    unit_test(test_y1),
//    unit_test(test_yn),
//    unit_test(test_dy0dx),
//    unit_test(test_dy1dx),
//    unit_test(test_dyndx)
  };
  return run_tests(tests);
}
