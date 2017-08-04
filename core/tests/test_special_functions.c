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
#include "core/polymec.h"
#include "core/special_functions.h"

// If you're messing with Bessel functions, you probably want double 
// precision.
#if POLYMEC_HAVE_DOUBLE_PRECISION

static void test_tgamma(void** state)
{
  assert_true(reals_nearly_equal(tgamma(-1.5), 2.363271801207, 1e-8));
  assert_true(reals_nearly_equal(tgamma(-0.5), -3.544907701811, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 0.5), 1.772453850906, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 1.0), 1.0, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 1.5), 0.886226925453, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 2.0), 1.0, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 2.5), 1.329340388179, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 3.0), 2.0, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 3.5), 3.323350970448, 1e-8));
  assert_true(reals_nearly_equal(tgamma( 4.0), 6.0, 1e-8));
}

static void test_bessel_find_jn_roots(void** state)
{
  static real_t jn_roots[6][5] = 
    {{2.4048, 5.5201, 8.6537, 11.7915, 14.9309},    // J0
     {3.8317, 7.0156, 10.1735, 13.3237, 16.4706},   // J1
     {5.1356, 8.4172, 11.6198, 14.7960, 17.9598},   // J2
     {6.3802, 9.7610, 13.0152, 16.2235, 19.4094},   // J3
     {7.5883, 11.0647, 14.3725, 17.6160, 20.8269},  // J4
     {8.7715, 12.3386, 15.7002, 18.9801, 22.2178}}; // J5
   
  for (int n = 0; n <= 5; ++n)
  {
    real_t roots[5];
    bessel_find_jn_roots(n, 5, roots);
    for (int k = 0; k < 5; ++k)
    {
      assert_true(ABS(roots[k] - jn_roots[n][k])/jn_roots[n][k] < 1e-4);
    }
  }
}

static void test_bessel_jn(void** state, int n)
{
  int num_roots = 10;
  real_t roots[num_roots];
  bessel_find_jn_roots(n, num_roots, roots);
  real_t sign = -1.0;
  for (int i = 0; i < num_roots; ++i)
  {
    assert_true(reals_nearly_equal(bessel_jn(n, roots[i]), 0.0, sqrt(REAL_EPSILON)));
    if (i < (num_roots-1))
    {
      assert_true(sign * bessel_jn(n, 0.5*(roots[i+1]+roots[i])) > 0.0);
      sign *= -1.0;
    }
  }
}

static void test_bessel_j(void** state)
{
  assert_true(reals_equal(bessel_j0(0.0), 1.0));
  assert_true(reals_equal(bessel_j1(0.0), 0.0));
  for (int n = 0; n < 20; ++n)
    test_bessel_jn(state, n);
}

static void test_bessel_djndx(void** state, int n)
{
  int num_roots = 10;
  real_t roots[num_roots];
  bessel_find_jn_roots(n, num_roots, roots);
  real_t sign = -1.0;
  for (int i = 0; i < num_roots; ++i)
  {
    assert_true(sign * bessel_djndx(n, roots[i]-0.01) > 0.0);
    sign *= -1.0;
  }
}

static void test_bessel_djdx(void** state)
{
  for (int n = 0; n < 20; ++n)
    test_bessel_djndx(state, n);
}

static void test_bessel_find_yn_roots(void** state)
{
  // These reference roots were taken from scipy.special.yn_zeros.
  static real_t yn_roots[6][5] = 
    {{0.8936,   3.9577,   7.0861,  10.2223,  13.3611},  // Y0
     {2.1971,   5.4297,   8.5960,  11.7492,  14.8974},  // Y1
     {3.3842,   6.7938,  10.0235,  13.2100,  16.3790},  // Y2
     {4.5270,   8.0976,  11.3965,  14.6231,  17.8185},  // Y3
     {5.6451,   9.3616,  12.7301,  15.9996,  19.2244},  // Y4
     {6.7472,  10.5972,  14.0338,  17.3471,  20.6029}}; // Y5

  for (int n = 0; n <= 5; ++n)
  {
    real_t roots[5];
    bessel_find_yn_roots(n, 5, roots);
    for (int k = 0; k < 5; ++k)
    {
      assert_true(ABS(roots[k] - yn_roots[n][k])/yn_roots[n][k] < 1e-4);
    }
  }
}

static void test_bessel_yn(void** state, int n)
{
  int num_roots = 10;
  real_t roots[num_roots];
  bessel_find_yn_roots(n, num_roots, roots);
  real_t sign = 1.0;
  for (int i = 0; i < num_roots; ++i)
  {
    assert_true(reals_nearly_equal(bessel_yn(n, roots[i]), 0.0, sqrt(REAL_EPSILON)));
    if (i < (num_roots-1))
    {
      assert_true(sign * bessel_yn(n, 0.5*(roots[i+1]+roots[i])) > 0.0);
      sign *= -1.0;
    }
  }
}

static void test_bessel_y(void** state)
{
  for (int n = 0; n < 20; ++n)
    test_bessel_yn(state, n);
}

static void test_bessel_dyndx(void** state, int n)
{
  int num_roots = 10;
  real_t roots[num_roots];
  bessel_find_yn_roots(n, num_roots, roots);
  real_t sign = 1.0;
  for (int i = 0; i < num_roots; ++i)
  {
    assert_true(sign * bessel_dyndx(n, roots[i]-0.01) > 0.0);
    sign *= -1.0;
  }
}

static void test_bessel_dydx(void** state)
{
  for (int n = 0; n < 20; ++n)
    test_bessel_dyndx(state, n);
}

static void test_chebyshev_tn(void** state)
{
  // Test the first 11 Chebyshev polynomials against those listed on 
  // Wikipedia.
  real_t x = 1.5; // an interesting point.
  assert_true(reals_nearly_equal(chebyshev_tn(0, x),  1.0, 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(1, x),  x, 1e-14));
printf("t2(1.5) == %f vs %f\n", chebyshev_tn(2, x), 2.0*x*x-1.0);
  assert_true(reals_nearly_equal(chebyshev_tn(2, x), (2.0*x*x-1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(3, x), (4.0*pow(x, 3)-3.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(4, x), (8.0*pow(x, 4)-8.0*x*x+1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(5, x), (16.0*pow(x, 5)-20.0*pow(x,3)+5.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(6, x), (32.0*pow(x, 6)-48.0*pow(x,4)+18.0*x*x-1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(7, x), (64.0*pow(x, 7)-112.0*pow(x, 5)+56.0*pow(x,3)-7.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(8, x), (128.0*pow(x, 8)-256.0*pow(x, 6)+160.0*pow(x,4)-32.0*pow(x,2)+1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(9, x), (256.0*pow(x, 8)-576.0*pow(x, 6)+432.0*pow(x,4)-120.0*pow(x,2)+9.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(10, x), (512.0*pow(x, 9)-1280.0*pow(x, 7)+1120.0*pow(x,5)-400.0*pow(x,3)+50.0*x*x-1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_tn(11, x), (1024.0*pow(x, 9)-2816.0*pow(x, 7)+2816.0*pow(x,5)-1232.0*pow(x,3)+220.0*x*x*x-11.0*x), 1e-14));
}

static void test_chebyshev_un(void** state)
{
  // Test the first 10 Chebyshev polynomials against those listed on 
  // Wikipedia.
  real_t x = 1.5; // an interesting point.
  assert_true(reals_nearly_equal(chebyshev_un(0, x), (1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(1, x), (2.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(2, x), (4.0*x*x-2.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(3, x), (8.0*pow(x, 3)-4.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(4, x), (16.0*pow(x, 4)-12.0*x*x+1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(5, x), (32.0*pow(x, 5)-32.0*pow(x,3)+6.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(6, x), (64.0*pow(x, 6)-80.0*pow(x,4)+24.0*x*x-1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(7, x), (128.0*pow(x, 7)-192.0*pow(x, 5)+80.0*pow(x,3)-8.0*x), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(8, x), (256.0*pow(x, 8)-448.0*pow(x, 6)+240.0*pow(x,4)-40.0*pow(x,2) + 1.0), 1e-14));
  assert_true(reals_nearly_equal(chebyshev_un(9, x), (512.0*pow(x, 9)-1024.0*pow(x, 7)+672.0*pow(x,5)-160.0*pow(x,3) + 10.0*x), 1e-14));
}

#endif // if double precision

static void test_hermite_hn(void** state)
{
  // Test the first 11 Hermite polynomials against those listed on 
  // Wikipedia.
  real_t x = 1.5; // an interesting point.
  assert_true(reals_nearly_equal(hermite_hn(0, x), (1.0), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(1, x), (2.0*x), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(2, x), (4.0*x*x-2.0), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(3, x), (8.0*pow(x, 3)-12.0*x), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(4, x), (16.0*pow(x, 4)-48.0*x*x+12.0), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(5, x), (32.0*pow(x, 5)-160.0*pow(x,3)+120.0*x), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(6, x), (64.0*pow(x, 6)-480.0*pow(x,4)+720.0*x*x-120.0), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(7, x), (128.0*pow(x, 7)-1344.0*pow(x, 5)+3360.0*pow(x,3)-1680.0*x), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(8, x), (256.0*pow(x, 8)-3584.0*pow(x, 6)+13440.0*pow(x,4)-13440.0*pow(x,2) + 1680.0), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(9, x), (512.0*pow(x, 9)-9216.0*pow(x, 7)+48384.0*pow(x,5)-80640.0*pow(x,3) + 30240.0*x), 1e-14));
  assert_true(reals_nearly_equal(hermite_hn(10, x), (1024.0*pow(x, 10)-23040.0*pow(x, 8)+161280.0*pow(x,6)-403200.0*pow(x,4) + 302400.0*x*x - 30240.0), 1e-14));
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
#if POLYMEC_HAVE_DOUBLE_PRECISION
    cmocka_unit_test(test_tgamma),
    cmocka_unit_test(test_bessel_find_jn_roots),
    cmocka_unit_test(test_bessel_j),
    cmocka_unit_test(test_bessel_djdx),
    cmocka_unit_test(test_bessel_find_yn_roots),
    cmocka_unit_test(test_bessel_y),
    cmocka_unit_test(test_bessel_dydx),
    cmocka_unit_test(test_chebyshev_tn),
    cmocka_unit_test(test_chebyshev_un),
#endif
    cmocka_unit_test(test_hermite_hn)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
