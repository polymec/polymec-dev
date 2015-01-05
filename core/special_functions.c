// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/special_functions.h"

// Much of the source code in this file originated from Fortran 77 subroutines 
// accompanying "Computation of Special Functions" by Shanjie Zhang and 
// Jianming Jin (Copyright 1996 by John Wiley & Sons, Inc). The Fortran code 
// was downloaded from http://jin.ece.illinois.edu/routines/routines.html.
// The C transcription is subject to Polymec's license, and the original 
// copyright should be mentioned in any distribution.

double gamma(double x)
{
  ASSERT((x >= 0.0) || (x != 1.0 * (int)x));
  ASSERT(x <= 171.6);

  // If x is an integer, we can do the easy thing.
  if (x == 1.0*(int)x)
  {
    double ga = 1.0;
    for (int i = 2; i < (int)x; ++i)
      ga *= i;
    return ga;
  }

  double z, r = 1.0;
  int m;
  if (fabs(x) > 1.0)
  {
    z = fabs(x);
    m = (int)z;
    for (int i = 1; i <= m; ++i)
      r *= (z - i);
    z -= 1.0*m;
  }
  else
    z = x;

  static double g[26] = { 1.0, 0.5772156649015329, 
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
  double gr = g[25];
  for (int i = 0; i < 25; ++i)
     gr = gr*z + g[24-i];
  double ga = 1.0 / (gr*z);
  if (fabs(x) > 1.0)
  {
    ga *= r;
    if (x < 0.0) 
      ga = -M_PI / (x * ga * sin(M_PI*x));
  }
  return ga;
}

// These data are used in the computation of J0/Y0 and J1/Y1.
static const double jy_rp2 = 0.63661977236758;
static double jy0_data_a[12] = 
  {-0.7031250000000000e-01, 0.1121520996093750e+00,
   -0.5725014209747314e+00, 0.6074042001273483e+01,
   -0.1100171402692467e+03, 0.3038090510922384e+04,
   -0.1188384262567832e+06, 0.6252951493434797e+07,
   -0.4259392165047669e+09, 0.3646840080706556e+11,
   -0.3833534661393944e+13, 0.4854014686852901e+15};
static double jy0_data_b[12] = 
  { 0.7324218750000000e-01,-0.2271080017089844e+00,
    0.1727727502584457e+01,-0.2438052969955606e+02,
    0.5513358961220206e+03,-0.1825775547429318e+05,
    0.8328593040162893e+06,-0.5006958953198893e+08,
    0.3836255180230433e+10,-0.3649010818849833e+12,
    0.4218971570284096e+14,-0.5827244631566907e+16};
static double jy1_data_a[12] = 
  { 0.1171875000000000e+00,-0.1441955566406250e+00,
    0.6765925884246826e+00,-0.6883914268109947e+01,
    0.1215978918765359e+03,-0.3302272294480852e+04,
    0.1276412726461746e+06,-0.6656367718817688e+07,
    0.4502786003050393e+09,-0.3833857520742790e+11,
    0.4011838599133198e+13,-0.5060568503314727e+15};
static double jy1_data_b[12] = 
  {-0.1025390625000000e+00, 0.2775764465332031e+00,
   -0.1993531733751297e+01, 0.2724882731126854e+02,
   -0.6038440767050702e+03, 0.1971837591223663e+05,
   -0.8902978767070678e+06, 0.5310411010968522e+08,
   -0.4043620325107754e+10, 0.3827011346598605e+12,
   -0.4406481417852278e+14, 0.6065091351222699e+16};

double j0(double x)
{
  if (x == 0.0) 
    return 1.0;
  
  double x2 = x*x;
  double val;
  if (x < 12.0)
  {
    val = 1.0;
    double r = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      r = -0.25*r*x2/(k*k);
      val += r;
      if (fabs(r) < 1e-15*fabs(val)) break;
    }
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    double t1 = x - 0.25*M_PI;
    double p0 = 1.0;
    double q0 = -0.125/x;
    for (int k = 1; k <= k0; ++k)
    {
      p0 += jy0_data_a[k-1] * pow(x, -2.0*k);
      q0 += jy0_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    double cu = sqrt(jy_rp2/x);
    val = cu * (p0*cos(t1) - q0*sin(t1));
  }
  return val;
}

double j1(double x)
{
  if (x == 0.0) 
    return 0.0;
  
  double x2 = x*x;
  double val;
  if (x < 12.0)
  {
    val = 1.0;
    double r = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      r = -0.25*r*x2/(k*(k+1));
      val += r;
      if (fabs(r) < 1e-15*fabs(val)) break;
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
    double t2 = x - 0.75*M_PI;
    double p1 = 1.0;
    double q1 = 0.375/x;
    for (int k = 1; k <= k0; ++k)
    {
      p1 += jy1_data_a[k-1] * pow(x, -2.0*k);
      q1 += jy1_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    double cu = sqrt(jy_rp2/x);
    val = cu * (p1*cos(t2) - q1*sin(t2));
  }
  return val;
}

// These functions are used for the calculation of Jn and Yn.
static inline double envj(int n, double x)
{
  return 0.5 * log10(6.28*n) - n * log10(1.36*x/n);
}

static double msta1(double x, double mp)
{
  double a0 = fabs(x);
  int n0 = (int)(1.1*a0) + 1;
  double f0 = envj(n0, a0) - mp;
  int n1 = n0 + 5;
  double f1 = envj(n1, a0) - mp;
  double nn;
  for (int it = 1; it <= 20; ++it)
  {
    nn = n1 - 1.0*(n1-n0)/(1.0 - f0/f1);
    double f = envj(nn, a0) - mp;
    if (fabs(nn-n1) < 1.0) break;
    n0 = n1;
    f0 = f1;
    n1 = nn;
    f1 = f;
  }
  return nn;
}

static double msta2(double x, int n, double mp)
{
  double a0 = fabs(x);
  double hmp = 0.5 * mp;
  double ejn = envj(n, a0);
  double obj;
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
  double f0 = envj(n0, a0) - obj;
  int n1 = n0 + 5;
  double f1 = envj(n1, a0) - obj;
  double nn;
  for (int it = 1; it <= 20; ++it)
  {
    nn = 1.0*n1 - 1.0*(n1 - n0)/(1.0 - f0/f1);
    double f = envj(nn, a0) - obj;
    if (fabs(nn-n1) < 1.0) break;
    n0 = n1;
    f0 = f1;
    n1 = nn;
    f1 = f;
  }
  return nn + 10;
}

double jn(int n, double x)
{
  if (n == 0)
    return j0(x);
  else if (n == 1)
    return j1(x);

  if (x < 1e-100)
    return 0.0;

  double val, j0_val = j0(x), j1_val = j1(x);
  double jl = j0_val, jm = j1_val;
  if (n < (int)(0.9*x))
  {
    for (int k = 2; k <= n; ++k)
    {
      double jk = 2.0 * jm * (k-1)/x - jl;
      val = jk;
      jl = jm;
      jm = jk;
    }
  }
  else
  {
    // Determine the order from which to recurse backwards.
    int nm = n;
    int m = msta1(x, 200.0);
    if (m < n)
      nm = m;
    else
      m = msta2(x, n, 15.0);

    double f, f2 = 0.0, f1 = 99.0;
    for (int k = m; k >= 0; --k)
    {
      f = 2.0*f1*(k+1)/x - f2;
      if (k == n) 
        val = f;
      f2 = f1;
      f1 = f;
    }
    double cs;
    if (fabs(j0_val) > fabs(j1_val))
      cs = j0_val / f;
    else
      cs = j1_val / f2;
    val *= cs;
  }

  return val;
}

double dj0dx(double x)
{
  return -j1(x);
}

double dj1dx(double x)
{
  return j0(x) - j1(x)/(x + 1e-15);
}

double djndx(int n, double x)
{
  if (n == 0)
    return dj0dx(x);
  else if (n == 1)
    return dj1dx(x);
  else
    return jn(n-1, x) - (1.0*(n+1)/x) * jn(n, x);
}

void find_jn_roots(int n, int num_roots, double* roots)
{
  // Find a reasonable starting point.
  double x;
  if (n <= 20) 
    x = 2.82141 + 1.15859 * n;
  else
  {
    double n_third = pow(n, 0.33333);
    x = 1.0*n + 1.85576 * n_third + 1.03315/n_third;
  }

  // Find the roots by Newton iteration.
  double x0;
  for (int i = 0; i < num_roots; ++i)
  {
    do
    {
      x0 = x;
      double jn_val = jn(n, x), djndx_val = djndx(n, x);
      x -= jn_val/djndx_val;
    }
    while (fabs(x-x0) > 1e-9);
    roots[i] = x;
    x += M_PI + (0.0972 + 0.0679*n - 0.000354*n*n)/(i+1);
  }
}

double y0(double x)
{
  if (x == 0.0) 
    return -1.0e300;
  
  double x2 = x*x;
  double val;
  if (x < 12.0)
  {
    double ec = log(0.5*x) + 0.5772156649015329;
    double cs0 = 0.0, w0 = 0.0, r0 = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      w0 += 1.0/k;
      r0 *= -0.25*x2 / (k*k);
      double r = r0*w0;
      cs0 += r;
      if (fabs(r) < 1e-15*fabs(cs0)) break;
    }
    double j0_val = j0(x);
    val = jy_rp2 * (ec*j0_val - cs0);
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    double t1 = x - 0.25*M_PI;
    double p0 = 1.0;
    double q0 = -0.125/x;
    for (int k = 1; k <= k0; ++k)
    {
      p0 += jy0_data_a[k-1] * pow(x, -2.0*k);
      q0 += jy0_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    double cu = sqrt(jy_rp2/x);
    val = cu * (p0*sin(t1) + q0*cos(t1));
  }
  return val;
}

double y1(double x)
{
  if (x == 0.0) 
    return -1.0e300;
  
  double x2 = x*x;
  double val;
  if (x < 12.0)
  {
    double ec = log(0.5*x) + 0.5772156649015329;
    double cs1 = 1.0, w1 = 0.0, r1 = 1.0;
    for (int k = 1; k <= 30; ++k)
    {
      w1 += 1.0/k;
      r1 *= -0.25*x2 / (k*(k+1));
      double r = r1 * (2.0*w1 + 1.0/(k+1));
      cs1 += r;
      if (fabs(r) < 1e-15*fabs(cs1)) break;
    }
    double j1_val = j1(x);
    val = jy_rp2 * (ec*j1_val - 1.0/x - 0.25*x*cs1);
  }
  else
  {
    int k0 = 12;
    if (x >= 50.0) 
      k0 = 10;
    else if (x >= 35.0)
      k0 = 8;
    double t2 = x - 0.75*M_PI;
    double p1 = 1.0;
    double q1 = 0.375/x;
    for (int k = 1; k <= k0; ++k)
    {
      p1 += jy1_data_a[k-1] * pow(x, -2.0*k);
      q1 += jy1_data_b[k-1] * pow(x, -2.0*k-1.0);
    }
    double cu = sqrt(jy_rp2/x);
    val = cu * (p1*sin(t2) + q1*cos(t2));
  }
  return val;
}

double yn(int n, double x)
{
  if (n == 0)
    return y0(x);
  else if (n == 1)
    return y1(x);

  if (x < 1e-100)
    return -1e300;

  double val, y0_val = y0(x), y1_val = y1(x);
  double f0 = y0_val, f1 = y1_val;
  for (int k = 2; k <= n; ++k)
  {
    double f = 2.0*f1*(k-1)/x - f0;
    if (k == n)
      val = f;
    f0 = f1;
    f1 = f;
  }
  return val;
}

double dy0dx(double x)
{
  return -y1(x);
}

double dy1dx(double x)
{
  return y0(x) - y1(x)/(x + 1e-15);
}

double dyndx(int n, double x)
{
  if (n == 0)
    return dy0dx(x);
  else if (n == 1)
    return dy1dx(x);
  else
    return yn(n-1, x) - (1.0*(n+1)/x) * yn(n, x);
}

void find_yn_roots(int n, int num_roots, double* roots)
{
  // Find a reasonable starting point.
  double x;
  if (n <= 20) 
    x = 1.19477 + 1.08933 * n;
  else
  {
    double n_third = pow(n, 0.33333);
    x = n + 0.93158 * n_third + 0.26035 / n_third;
  }

  // Find the roots by Newton iteration.
  double x0;
  for (int i = 0; i < num_roots; ++i)
  {
    do
    {
      x0 = x;
      double yn_val = yn(n, x), dyndx_val = dyndx(n, x);
      x -= yn_val/dyndx_val;
    }
    while (fabs(x-x0) > 1e-9);
    roots[i] = x;
    x += M_PI + (0.312 + 0.0852*n - 0.000403*n*n)/(i+1);
  }
}

double i0(double x)
{
  POLYMEC_NOT_IMPLEMENTED
}

double i1(double x)
{
  POLYMEC_NOT_IMPLEMENTED
}

double in(int n, double x)
{
  POLYMEC_NOT_IMPLEMENTED
}

double k0(double x)
{
  POLYMEC_NOT_IMPLEMENTED
}

double k1(double x)
{
  POLYMEC_NOT_IMPLEMENTED
}

double kn(int n, double x)
{
  POLYMEC_NOT_IMPLEMENTED
}

#ifndef __cplusplus
#include <complex.h>

double complex cj0(double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex cj1(double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex cjn(int n, double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex cy0(double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex cy1(double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex cyn(int n, double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex ci0(double complex z)
{
  double a0 = cabs(z);
  double complex i0 = CMPLX(1.0, 0.0);
  if (a0 == 0.0)
    return i0;
  double complex z1 = (creal(z) < 0.0) ? -z : z;
  double complex z2 = z * z;
  if (a0 <= 18.0)
  {
    double complex cr = CMPLX(1.0, 0.0);
    for (int i = 1; i <= 50; ++i)
    {
      cr = 0.25 * cr * z2 / (i * i);
      i0 += cr;
      if (cabs(cr/i0) < 1e-15) break;
    }
  }
  else
  {
    static double a[12] = { 0.125,              7.03125e-2,
                            7.32421875e-2,      1.1215209960938e-1,
                            2.2710800170898e-1, 5.7250142097473e-1,
                            1.7277275025845,    6.0740420012735,
                            2.4380529699556e1,  1.1001714026925e2,
                            5.5133589612202e2,  3.0380905109224e3};
    double k0 = 12.0;
    if (a0 >= 50.0)
      k0 = 7.0;
    else if (a0 >= 35.0) 
      k0 = 9.0;
    double complex ca = cexp(z1) / csqrt(2.0 * M_PI * z1);
    double complex zr = 1.0 / z1;
    for (int i = 1; i <= (int)k0; ++i)
      i0 += a[i-1] * cpow(zr, i);
    i0 *= ca;
  }
  return i0;
}

double complex ci1(double complex z)
{
  double a0 = cabs(z);
  if (a0 == 0.0)
    return CMPLX(0.0, 0.0);
  double complex z1 = (creal(z) < 0.0) ? -z : z;
  double complex z2 = z * z;
  double complex i1 = CMPLX(1.0, 0.0);
  if (a0 <= 18.0)
  {
    double complex cr = CMPLX(1.0, 0.0);
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
    static double b[12] = {-0.375,              -1.171875e-1,
                           -1.025390625e-1,     -1.4419555664063e-1,
                           -2.7757644653320e-1, -6.7659258842468e-1,
                           -1.9935317337513,    -6.8839142681099,
                           -2.7248827311269e1,  -1.2159789187654e2,
                           -6.0384407670507e2,  -3.3022722944809e3};
    double k0 = 12.0;
    if (a0 >= 50.0)
      k0 = 7.0;
    else if (a0 >= 35.0) 
      k0 = 9.0;
    double complex ca = cexp(z1) / csqrt(2.0 * M_PI * z1);
    double complex zr = 1.0 / z1;
    for (int i = 1; i <= (int)k0; ++i)
      i1 += b[i-1] * cpow(zr, i);
    i1 *= ca;
  }
  if (creal(z) < 0.0)
    i1 = -i1;
  return i1;
}

double complex cin(int n, double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

double complex ck0(double complex z)
{
  double a0 = cabs(z);
  if (a0 == 0.0)
    return CMPLX(1e300, 0.0);
  double complex z1 = (creal(z) < 0.0) ? -z : z;
  double complex z2 = z * z;
  double complex i0 = ci0(z);
  double complex k0;
  if (a0 <= 9.0)
  {
    double complex cs = CMPLX(0.0, 0.0);
    double complex ct = -clog(0.5*z1) - 0.5772156649015329;
    double w0 = 0.0;
    double complex cr = CMPLX(1.0, 0.0);
    double complex cw = cs;
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
    static double a1[12] = {0.125,             0.2109375,
                            1.0986328125,      1.1775970458984e1,
                            2.1461706161499e2, 5.9511522710323e3,
                            2.3347645606175e5, 1.2312234987631e7,
                            8.401390346421e8,  7.2031420482627e10};
    double complex cb = 0.5/z1;
    double complex zr2 =1.0/z2;
    k0 = CMPLX(1.0, 0.0);
    for (int i = 1; i <= 10; ++i)
      k0 += a1[i-1] * cpow(zr2, i);
    k0 = cb * k0 / i0;
  }
  if (creal(z) < 0.0)
  {
    if (cimag(z) < 0.0)
      k0 += I * M_PI * i0;
    else if (cimag(z) > 0.0)
      k0 -= I * M_PI * i0;
  }
  return k0;
}

double complex ck1(double complex z)
{
  double a0 = cabs(z);
  if (a0 == 0.0)
    return CMPLX(1e300, 0.0);
  double complex z1 = (creal(z) < 0.0) ? -z : z;
  double complex i0 = ci0(z);
  double complex i1 = ci1(z);
  double complex k0 = ck0(z);
  double complex k1 = (1.0/z1 - i1*k0) / i0;
  if (creal(z) < 0.0)
  {
    if (cimag(z) < 0.0)
      k1 = -k1 + I * M_PI * i1;
    else if (cimag(z) > 0.0)
      k1 = -k1 - I * M_PI * i1;
  }
  return k1;
}

double complex ckn(int n, double complex z)
{
  POLYMEC_NOT_IMPLEMENTED
}

#endif

