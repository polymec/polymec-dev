// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/special_functions.h"

#ifdef _Imaginary_I
#define I _Imaginary_I
#else
#define I _Complex_I
#endif

// FIXME: Remove these pragmas when all special functions are implemented.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsuggest-attribute=noreturn"

// Much of the source code in this file originated from Fortran 77 subroutines 
// accompanying "Computation of Special Functions" by Shanjie Zhang and 
// Jianming Jin (Copyright 1996 by John Wiley & Sons, Inc). The Fortran code 
// was downloaded from http://jin.ece.illinois.edu/routines/routines.html.
// The C transcription is subject to Polymec's license, and the original 
// copyright should be mentioned in any distribution.

real_t gamma(real_t x)
{
  ASSERT((x >= 0.0) || reals_equal(x, 1.0*(int)x));
  ASSERT(x <= 171.6);

  // If x is an integer, we can do the easy thing.
  if (reals_equal(x, 1.0*(int)x))
  {
    real_t ga = 1.0;
    for (int i = 2; i < (int)x; ++i)
      ga *= i;
    return ga;
  }

  real_t z, r = 1.0;
  int m;
  if (ABS(x) > 1.0)
  {
    z = ABS(x);
    m = (int)z;
    for (int i = 1; i <= m; ++i)
      r *= (z - i);
    z -= 1.0*m;
  }
  else
    z = x;

  static real_t g[26] = { 1.0, 0.5772156649015329, 
                         -0.6558780715202538,  -0.420026350340952e-1,
                          0.1665386113822915,  -0.421977345555443e-1, 
                         -0.9621971527877e-2,   0.7218943246663e-2,
                         -0.11651675918591e-2, -0.2152416741149e-3,
                          0.1280502823882e-3,  -0.201348547807e-4,
                         -0.12504934821e-5, 0.1133027232e-5,
                         -0.2056338417e-6, 0.6116095e-8,
                          0.50020075e-8, -0.11812746e-8,
                          0.1043427e-9, 0.77823e-11,
                         -0.36968e-11, 0.51e-12,
                         -0.206e-13, -0.54e-14,
                          0.14e-14, 0.1e-14};
  real_t gr = g[25];
  for (int i = 0; i < 25; ++i)
     gr = gr*z + g[24-i];
  real_t ga = 1.0 / (gr*z);
  if (ABS(x) > 1.0)
  {
    ga *= r;
    if (x < 0.0) 
      ga = -M_PI / (x * ga * sin(M_PI*x));
  }
  return ga;
}

// These data are used in the computation of J0/Y0 and J1/Y1.
static const real_t jy_rp2 = 0.63661977236758;
static real_t jy0_data_a[12] = 
  {-0.7031250000000000e-01, 0.1121520996093750e+00,
   -0.5725014209747314e+00, 0.6074042001273483e+01,
   -0.1100171402692467e+03, 0.3038090510922384e+04,
   -0.1188384262567832e+06, 0.6252951493434797e+07,
   -0.4259392165047669e+09, 0.3646840080706556e+11,
   -0.3833534661393944e+13, 0.4854014686852901e+15};
static real_t jy0_data_b[12] = 
  { 0.7324218750000000e-01,-0.2271080017089844e+00,
    0.1727727502584457e+01,-0.2438052969955606e+02,
    0.5513358961220206e+03,-0.1825775547429318e+05,
    0.8328593040162893e+06,-0.5006958953198893e+08,
    0.3836255180230433e+10,-0.3649010818849833e+12,
    0.4218971570284096e+14,-0.5827244631566907e+16};
static real_t jy1_data_a[12] = 
  { 0.1171875000000000e+00,-0.1441955566406250e+00,
    0.6765925884246826e+00,-0.6883914268109947e+01,
    0.1215978918765359e+03,-0.3302272294480852e+04,
    0.1276412726461746e+06,-0.6656367718817688e+07,
    0.4502786003050393e+09,-0.3833857520742790e+11,
    0.4011838599133198e+13,-0.5060568503314727e+15};
static real_t jy1_data_b[12] = 
  {-0.1025390625000000e+00, 0.2775764465332031e+00,
   -0.1993531733751297e+01, 0.2724882731126854e+02,
   -0.6038440767050702e+03, 0.1971837591223663e+05,
   -0.8902978767070678e+06, 0.5310411010968522e+08,
   -0.4043620325107754e+10, 0.3827011346598605e+12,
   -0.4406481417852278e+14, 0.6065091351222699e+16};

real_t bessel_j0(real_t x)
{
  if (reals_equal((real_t)x, 0.0)) 
    return 1.0;
  
  real_t x2 = x*x;
  real_t val;
  if (x < 12.0)
  {
    val = 1.0;
    real_t r = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      r = -0.25*r*x2/(k*k);
      val += r;
      if (ABS(r) < 1e-15*ABS(val)) break;
    }
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    real_t t1 = x - 0.25*M_PI;
    real_t p0 = 1.0;
    real_t q0 = -0.125/x;
    for (int k = 1; k <= k0; ++k)
    {
      p0 += jy0_data_a[k-1] * pow(x, -2.0*k);
      q0 += jy0_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    real_t cu = sqrt(jy_rp2/x);
    val = cu * (p0*cos(t1) - q0*sin(t1));
  }
  return val;
}

real_t bessel_j1(real_t x)
{
  if (reals_equal((real_t)x, 0.0)) 
    return 0.0;
  
  real_t x2 = x*x;
  real_t val;
  if (x < 12.0)
  {
    val = 1.0;
    real_t r = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      r = -0.25*r*x2/(k*(k+1));
      val += r;
      if (ABS(r) < 1e-15*ABS(val)) break;
    }
    val *= 0.5*x;
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    real_t t2 = x - 0.75*M_PI;
    real_t p1 = 1.0;
    real_t q1 = 0.375/x;
    for (int k = 1; k <= k0; ++k)
    {
      p1 += jy1_data_a[k-1] * pow(x, -2.0*k);
      q1 += jy1_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    real_t cu = sqrt(jy_rp2/x);
    val = cu * (p1*cos(t2) - q1*sin(t2));
  }
  return val;
}

// These functions are used for the calculation of Jn and Yn.
static inline real_t envj(int n, real_t x)
{
  return 0.5 * log10(6.28*n) - n * log10(1.36*x/n);
}

static real_t msta1(real_t x, real_t mp)
{
  real_t a0 = ABS(x);
  int n0 = (int)(1.1*a0) + 1;
  real_t f0 = envj(n0, a0) - mp;
  int n1 = n0 + 5;
  real_t f1 = envj(n1, a0) - mp;
  real_t nn = 0.0;
  for (int it = 1; it <= 20; ++it)
  {
    nn = n1 - 1.0*(n1-n0)/(1.0 - f0/f1);
    real_t f = envj((int)nn, a0) - mp;
    if (ABS(nn-n1) < 1.0) break;
    n0 = n1;
    f0 = f1;
    n1 = (int)nn;
    f1 = f;
  }
  return nn;
}

static real_t msta2(real_t x, int n, real_t mp)
{
  real_t a0 = ABS(x);
  real_t hmp = 0.5 * mp;
  real_t ejn = envj(n, a0);
  real_t obj;
  int n0;
  if (ejn <= hmp)
  {
    obj = mp;
    n0 = (int)(1.1*a0);
  }
  else
  {
    obj = hmp + ejn;
    n0 = n;
  }
  real_t f0 = envj(n0, a0) - obj;
  int n1 = n0 + 5;
  real_t f1 = envj(n1, a0) - obj;
  real_t nn = 0.0;
  for (int it = 1; it <= 20; ++it)
  {
    nn = 1.0*n1 - 1.0*(n1 - n0)/(1.0 - f0/f1);
    real_t f = envj((int)nn, a0) - obj;
    if (ABS(nn-n1) < 1.0) break;
    n0 = n1;
    f0 = f1;
    n1 = (int)nn;
    f1 = f;
  }
  return nn + 10.0;
}

real_t bessel_jn(int n, real_t x)
{
  if (n == 0)
    return bessel_j0(x);
  else if (n == 1)
    return bessel_j1(x);

  if (x < 1e-100)
    return 0.0;

  real_t val = 0.0, j0_val = bessel_j0(x), j1_val = bessel_j1(x);
  real_t jl = j0_val, jm = j1_val;
  if (n < (int)(0.9*x))
  {
    for (int k = 2; k <= n; ++k)
    {
      real_t jk = 2.0 * jm * (k-1)/x - jl;
      val = jk;
      jl = jm;
      jm = jk;
    }
  }
  else
  {
    // Determine the order from which to recurse backwards.
    int m = (int)(msta1(x, 200.0));
    if (m >= n)
      m = (int)(msta2(x, n, 15.0));

    real_t f = 0.0, f2 = 0.0, f1 = 99.0;
    for (int k = m; k >= 0; --k)
    {
      f = 2.0*f1*(k+1)/x - f2;
      if (k == n) 
        val = f;
      f2 = f1;
      f1 = f;
    }
    real_t cs;
    if (ABS(j0_val) > ABS(j1_val))
      cs = j0_val / f;
    else
      cs = j1_val / f2;
    val *= cs;
  }

  return val;
}

real_t bessel_dj0dx(real_t x)
{
  return -bessel_j1(x);
}

real_t bessel_dj1dx(real_t x)
{
  return bessel_j0(x) - bessel_j1(x)/(x + 1e-15);
}

real_t bessel_djndx(int n, real_t x)
{
  if (n == 0)
    return bessel_dj0dx(x);
  else if (n == 1)
    return bessel_dj1dx(x);
  else
    return bessel_jn(n-1, x) - (1.0*(n+1)/x) * bessel_jn(n, x);
}

void bessel_find_jn_roots(int n, 
                          int num_roots, 
                          real_t* roots)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  static const real_t tolerance = 1e-9;
#else
  static const real_t tolerance = 1e-4;
#endif

  // Find a reasonable starting point.
  real_t x;
  if (n <= 20) 
    x = 2.82141 + 1.15859 * n;
  else
  {
    real_t n_third = pow(n, 0.33333);
    x = 1.0*n + 1.85576 * n_third + 1.03315/n_third;
  }

  // Find the roots by Newton iteration.
  real_t x0;
  for (int i = 0; i < num_roots; ++i)
  {
    do
    {
      x0 = x;
      real_t jn_val = bessel_jn(n, x), djndx_val = bessel_djndx(n, x);
      x -= jn_val/djndx_val;
    }
    while (ABS(x-x0) > tolerance);
    roots[i] = x;
    x += M_PI + (0.0972 + 0.0679*n - 0.000354*n*n)/(i+1);
  }
}

real_t bessel_y0(real_t x)
{
  if (reals_equal((real_t)x, 0.0)) 
    return -1.0e300;
  
  real_t x2 = x*x;
  real_t val;
  if (x < 12.0)
  {
    real_t ec = log(0.5*x) + 0.5772156649015329;
    real_t cs0 = 0.0, w0 = 0.0, r0 = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      w0 += 1.0/k;
      r0 *= -0.25*x2 / (k*k);
      real_t r = r0*w0;
      cs0 += r;
      if (ABS(r) < 1e-15*ABS(cs0)) break;
    }
    real_t j0_val = bessel_j0(x);
    val = jy_rp2 * (ec*j0_val - cs0);
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    real_t t1 = x - 0.25*M_PI;
    real_t p0 = 1.0;
    real_t q0 = -0.125/x;
    for (int k = 1; k <= k0; ++k)
    {
      p0 += jy0_data_a[k-1] * pow(x, -2.0*k);
      q0 += jy0_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    real_t cu = sqrt(jy_rp2/x);
    val = cu * (p0*sin(t1) + q0*cos(t1));
  }
  return val;
}

real_t bessel_y1(real_t x)
{
  if (reals_equal((real_t)x, 0.0)) 
    return -1.0e300;
  
  real_t x2 = x*x;
  real_t val;
  if (x < 12.0)
  {
    real_t ec = log(0.5*x) + 0.5772156649015329;
    real_t cs1 = 1.0, w1 = 0.0, r1 = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      w1 += 1.0/k;
      r1 *= -0.25*x2 / (k*(k+1));
      real_t r = r1 * (2.0*w1 + 1.0/(k+1));
      cs1 += r;
      if (ABS(r) < 1e-15*ABS(cs1)) break;
    }
    real_t j1_val = bessel_j1(x);
    val = jy_rp2 * (ec*j1_val - 1.0/x - 0.25*x*cs1);
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    real_t t2 = x - 0.75*M_PI;
    real_t p1 = 1.0;
    real_t q1 = 0.375/x;
    for (int k = 1; k <= k0; ++k)
    {
      p1 += jy1_data_a[k-1] * pow(x, -2.0*k);
      q1 += jy1_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    real_t cu = sqrt(jy_rp2/x);
    val = cu * (p1*sin(t2) + q1*cos(t2));
  }
  return val;
}

real_t bessel_yn(int n, real_t x)
{
  if (n == 0)
    return bessel_y0(x);
  else if (n == 1)
    return bessel_y1(x);

  if (x < 1e-100)
    return -1e300;

  real_t val = 0.0, y0_val = bessel_y0(x), y1_val = bessel_y1(x);
  real_t f0 = y0_val, f1 = y1_val;
  for (int k = 2; k <= n; ++k)
  {
    real_t f = 2.0*f1*(k-1)/x - f0;
    if (k == n)
      val = f;
    f0 = f1;
    f1 = f;
  }
  return val;
}

real_t bessel_dy0dx(real_t x)
{
  return -bessel_y1(x);
}

real_t bessel_dy1dx(real_t x)
{
  return bessel_y0(x) - bessel_y1(x)/(x + 1e-15);
}

real_t bessel_dyndx(int n, real_t x)
{
  if (n == 0)
    return bessel_dy0dx(x);
  else if (n == 1)
    return bessel_dy1dx(x);
  else
    return bessel_yn(n-1, x) - (1.0*(n+1)/x) * bessel_yn(n, x);
}

void bessel_find_yn_roots(int n, 
                          int num_roots, 
                          real_t* roots)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  static const real_t tolerance = 1e-9;
#else
  static const real_t tolerance = 1e-4;
#endif

  // Find a reasonable starting point.
  real_t x;
  if (n <= 20) 
    x = 1.19477 + 1.08933 * n;
  else
  {
    real_t n_third = pow(n, 0.33333);
    x = n + 0.93158 * n_third + 0.26035 / n_third;
  }

  // Find the roots by Newton iteration.
  real_t x0;
  for (int i = 0; i < num_roots; ++i)
  {
    do
    {
      x0 = x;
      real_t yn_val = bessel_yn(n, x), dyndx_val = bessel_dyndx(n, x);
      x -= yn_val/dyndx_val;
    }
    while (ABS(x-x0) > tolerance);
    roots[i] = x;
    x += M_PI + (0.312 + 0.0852*n - 0.000403*n*n)/(i+1);
  }
}

real_t bessel_i0(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_i1(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_in(int n, real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_di0dx(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_di1dx(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_dindx(int n, real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_k0(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_k1(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_kn(int n, real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_dk0dx(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_dk1dx(real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t bessel_dkndx(int n, real_t x)
{
  POLYMEC_NOT_IMPLEMENTED
}

real_t chebyshev_tn(int n, real_t x)
{
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return x;
  real_t A = 2.0, B = 0.0, C = 1.0, y0 = 1.0, y1 = x;
  real_t tn = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    real_t yn = (A*x*B)*y1-C*y0;
    tn = yn;
    y0 = y1;
    y1 = yn;
  }
  return tn;
}

real_t chebyshev_dtndx(int n, real_t x)
{
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return 1.0;
  real_t A = 2.0, B = 0.0, C = 1.0, y0 = 1.0, y1 = x;
  real_t dy0 = 0.0, dy1 = 1.0;
  real_t dtndx = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    real_t yn = (A*x*B)*y1-C*y0;
    real_t dyn = A*y1+(A*x+B)*dy1-C*dy0;
    dtndx = dyn;
    y0 = y1;
    y1 = yn;
    dy0 = dy1;
    dy1 = dyn;
  }
  return dtndx;
}

real_t chebyshev_un(int n, real_t x)
{
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return 2.0*x;
  real_t A = 2.0, B = 0.0, C = 1.0, y0 = 1.0, y1 = 2.0*x;
  real_t un = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    real_t yn = (A*x*B)*y1-C*y0;
    un = yn;
    y0 = y1;
    y1 = yn;
  }
  return un;
}

real_t chebyshev_dundx(int n, real_t x)
{
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return 2.0;
  real_t A = 2.0, B = 0.0, C = 1.0, y0 = 1.0, y1 = 2.0*x;
  real_t dy0 = 0.0, dy1 = 2.0;
  real_t dundx = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    real_t yn = (A*x*B)*y1-C*y0;
    real_t dyn = A*y1+(A*x+B)*dy1-C*dy0;
    dundx = dyn;
    y0 = y1;
    y1 = yn;
    dy0 = dy1;
    dy1 = dyn;
  }
  return dundx;
}

real_t laguerre_ln(int n, real_t x)
{
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return 1.0 - x;
  real_t A = 2.0, B = 0.0, C = 1.0, y0 = 1.0, y1 = 1.0 - x;
  real_t ln = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    A = -1.0/k;
    B = 2.0 + A;
    C = 1.0 + A;
    real_t yn = (A*x*B)*y1-C*y0;
    ln = yn;
    y0 = y1;
    y1 = yn;
  }
  return ln;
}

real_t laguerre_dlndx(int n, real_t x)
{
  if (n == 0)
    return 0.0;
  else if (n == 1)
    return -1.0;
  real_t A = 2.0, B = 0.0, C = 1.0, y0 = 1.0, y1 = 1.0 - x;
  real_t dy0 = 0.0, dy1 = -1.0;
  real_t dlndx = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    A = -1.0/k;
    B = 2.0 + A;
    C = 1.0 + A;
    real_t yn = (A*x*B)*y1-C*y0;
    real_t dyn = A*y1+(A*x+B)*dy1-C*dy0;
    dlndx = dyn;
    y0 = y1;
    y1 = yn;
    dy0 = dy1;
    dy1 = dyn;
  }
  return dlndx;
}

real_t hermite_hn(int n, real_t x)
{
  ASSERT(n >= 0);
  if (n == 0)
    return 1.0;
  else if (n == 1)
    return 2.0*x;
  real_t h0 = 1.0, h1 = 2.0*x;
  real_t hn = 0.0;
  for (int k = 2; k <= n; ++k)
  {
    hn = 2.0*(x*h1 - (k-1)*h0);
    h0 = h1;
    h1 = hn;
  }
  return hn;
}

real_t hermite_dhndx(int n, real_t x)
{
  ASSERT(n >= 0);
  if (n == 0)
    return 0.0;
  else 
    return 2.0 * n * hermite_hn(n-1, x);
}

#ifndef __cplusplus
#include <complex.h>

complex_t bessel_h0(real_t x)
{
  return CMPLX(bessel_j0(x), bessel_y0(x));
}

complex_t bessel_h1(real_t x)
{
  return CMPLX(bessel_j1(x), bessel_y1(x));
}

complex_t bessel_hn(int n, real_t x)
{
  return CMPLX(bessel_jn(n, x), bessel_yn(n, x));
}

complex_t bessel_dh0dx(real_t x)
{
  return CMPLX(bessel_dj0dx(x), bessel_dy0dx(x));
}

complex_t bessel_dh1dx(real_t x)
{
  return CMPLX(bessel_dj1dx(x), bessel_dy1dx(x));
}

complex_t bessel_dhndx(int n, real_t x)
{
  return CMPLX(bessel_djndx(n, x), bessel_dyndx(n, x));
}

complex_t bessel_cj0(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_cj1(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_cjn(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_cy0(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_cy1(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_cyn(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_ci0(complex_t z)
{
  real_t a0 = cabs(z);
  complex_t i0 = CMPLX(1.0, 0.0);
  if (reals_equal((real_t)a0, 0.0)) 
    return i0;
  complex_t z1 = (creal(z) < 0.0) ? -z : z;
  complex_t z2 = z * z;
  if (a0 <= 18.0)
  {
    complex_t cr = CMPLX(1.0, 0.0);
    for (int i = 1; i <= 50; ++i)
    {
      cr = 0.25 * cr * z2 / (i * i);
      i0 += cr;
      if (cabs(cr/i0) < 1e-15) break;
    }
  }
  else
  {
    static real_t a[12] = { 0.125,              7.03125e-2,
                            7.32421875e-2,      1.1215209960938e-1,
                            2.2710800170898e-1, 5.7250142097473e-1,
                            1.7277275025845,    6.0740420012735,
                            2.4380529699556e1,  1.1001714026925e2,
                            5.5133589612202e2,  3.0380905109224e3};
    real_t k0 = 12.0;
    if (a0 >= 50.0)
      k0 = 7.0;
    else if (a0 >= 35.0) 
      k0 = 9.0;
    complex_t ca = cexp(z1) / csqrt(2.0 * M_PI * z1);
    complex_t zr = 1.0 / z1;
    for (int i = 1; i <= (int)k0; ++i)
      i0 += a[i-1] * cpow(zr, i);
    i0 *= ca;
  }
  return i0;
}

complex_t bessel_ci1(complex_t z)
{
  real_t a0 = cabs(z);
  if (reals_equal((real_t)a0, 0.0)) 
    return CMPLX(0.0, 0.0);
  complex_t z1 = (creal(z) < 0.0) ? -z : z;
  complex_t z2 = z * z;
  complex_t i1 = CMPLX(1.0, 0.0);
  if (a0 <= 18.0)
  {
    complex_t cr = CMPLX(1.0, 0.0);
    for (int i = 1; i <= 50; ++i)
    {
      cr = 0.25 * cr * z2 / (i * (i+1));
      i1 += cr;
      if (cabs(cr/i1) < 1e-15) break;
    }
    i1 *= 0.5*z1;
  }
  else
  {
    static real_t b[12] = {-0.375,              -1.171875e-1,
                           -1.025390625e-1,     -1.4419555664063e-1,
                           -2.7757644653320e-1, -6.7659258842468e-1,
                           -1.9935317337513,    -6.8839142681099,
                           -2.7248827311269e1,  -1.2159789187654e2,
                           -6.0384407670507e2,  -3.3022722944809e3};
    real_t k0 = 12.0;
    if (a0 >= 50.0)
      k0 = 7.0;
    else if (a0 >= 35.0) 
      k0 = 9.0;
    complex_t ca = cexp(z1) / csqrt(2.0 * M_PI * z1);
    complex_t zr = 1.0 / z1;
    for (int i = 1; i <= (int)k0; ++i)
      i1 += b[i-1] * cpow(zr, i);
    i1 *= ca;
  }
  if (creal(z) < 0.0)
    i1 = -i1;
  return i1;
}

complex_t bessel_cin(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_ck0(complex_t z)
{
  real_t a0 = cabs(z);
  if (reals_equal((real_t)a0, 0.0)) 
    return CMPLX(1e300, 0.0);
  complex_t z1 = (creal(z) < 0.0) ? -z : z;
  complex_t z2 = z * z;
  complex_t i0 = bessel_ci0(z);
  complex_t k0;
  if (a0 <= 9.0)
  {
    complex_t cs = CMPLX(0.0, 0.0);
    complex_t ct = -clog(0.5*z1) - 0.5772156649015329;
    real_t w0 = 0.0;
    complex_t cr = CMPLX(1.0, 0.0);
    complex_t cw = cs;
    for (int i = 1; i <= 50; ++i)
    {
      w0 += 1.0/i;
      cr *= 0.25 * z2 / (i*i);
      cs += cr * (w0 +ct);
      if (cabs((cs - cw)/cs) < 1e-15) 
        break;
      cw = cs;
    }
    k0 = ct + cs;
  }
  else
  {
    static real_t a1[12] = {0.125,             0.2109375,
                            1.0986328125,      1.1775970458984e1,
                            2.1461706161499e2, 5.9511522710323e3,
                            2.3347645606175e5, 1.2312234987631e7,
                            8.401390346421e8,  7.2031420482627e10};
    complex_t cb = 0.5/z1;
    complex_t zr2 =1.0/z2;
    k0 = CMPLX(1.0, 0.0);
    for (int i = 1; i <= 10; ++i)
      k0 += a1[i-1] * cpow(zr2, i);
    k0 = cb * k0 / i0;
  }
  if (creal(z) < 0.0)
  {
    complex_t i = CMPLX(0.0, 1.0);
    if (cimag(z) < 0.0)
      k0 += i * M_PI * i0;
    else if (cimag(z) > 0.0)
      k0 -= i * M_PI * i0;
  }
  return k0;
}

complex_t bessel_ck1(complex_t z)
{
  real_t a0 = cabs(z);
  if (reals_equal((real_t)a0, 0.0)) 
    return CMPLX(1e300, 0.0);
  complex_t z1 = (creal(z) < 0.0) ? -z : z;
  complex_t i0 = bessel_ci0(z);
  complex_t i1 = bessel_ci1(z);
  complex_t k0 = bessel_ck0(z);
  complex_t k1 = (1.0/z1 - i1*k0) / i0;
  if (creal(z) < 0.0)
  {
    complex_t i = CMPLX(0.0, 1.0);
    if (cimag(z) < 0.0)
      k1 = -k1 + i * M_PI * i1;
    else if (cimag(z) > 0.0)
      k1 = -k1 - i * M_PI * i1;
  }
  return k1;
}

complex_t bessel_ckn(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

// FIXME: Remove this when everything is implemented.
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

#endif

