// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/special_functions.h"

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

