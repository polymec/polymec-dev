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

real_t factorial(int n)
{
  ASSERT(n >= 0);

  static real_t fact[33] = {1.0, 1.0};
  static int top = 1;
  if (n > 32)
    return tgamma(1.0*(n+1));
  while (top < n)
  {
    fact[top+1] = top * fact[top];
    ++top;
  }
  return fact[n];
}

real_t binomial_coeff(int n, int k)
{
  ASSERT(k >= 0);
  ASSERT(k <= n);
  return floor(0.5+exp(lgamma(1.0*(n+1)) - lgamma(1.0*(k+1)) - lgamma(1.0*(n-k))));
}

real_t legendre_pn(int n, real_t x)
{
  real_t sum = 0.0;
  for (int k = 0; k <= n; ++k)
  {
    real_t n_choose_k = binomial_coeff(n, k);
    sum += pow(-1, k) * n_choose_k * n_choose_k * pow(0.5*(x+1.0), n-k) * pow(0.5*(x-1.0), k);
  }
  return sum;
}

real_t legendre_pml(int m, int l, real_t x)
{
  ASSERT(m >= 0);
  ASSERT(m <= l);
  ASSERT(x <= 1.0);

  real_t pmm = 1.0;
  if (m > 0)
  {
    real_t sq = (1.0-x) * (1.0+x);
    real_t fact = 1.0;
    for (int i = 1; i <= m; ++i)
    {
      pmm *= -fact * sq;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else
  {
    real_t pmmp = x * (2 * m + 1) * pmm;
    if (l == (m + 1))
      return pmmp;
    else
    {
      real_t pll = 0.0;
      for (int ll = m + 2; ll <= l; ++ll)
      {
        pll = (x * (2*ll-1) * pmmp - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp;
        pmmp = pll;
      }
      return pll;
    }
  }
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

complex_t bessel_jv(complex_t v, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dj0dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dj1dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_djndz(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_djvdz(complex_t v, complex_t z)
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

complex_t bessel_yv(complex_t v, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dy0dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dy1dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dyndz(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dyvdz(complex_t v, complex_t z)
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

complex_t bessel_iv(complex_t v, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_di0dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_di1dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

// Returns the value of the first derivative of In, the modified Bessel 
// function of the first kind of integer order n, at the given complex value z.
complex_t bessel_dindz(int n, complex_t z);

// Returns the value of the first derivative of Iv, the modified Bessel 
// function of the first kind of complex order nu (written v), at the given 
// complex value z.
complex_t bessel_divdz(complex_t v, complex_t z);
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

complex_t bessel_kv(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dk0dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dk1dz(complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dkndz(int n, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t bessel_dkvdz(complex_t v, complex_t z)
{
  POLYMEC_NOT_IMPLEMENTED
}

complex_t sph_harmonic_ylm(int l, int m, real_t theta, real_t phi)
{
  return pow(-1, m) * 
         sqrt((2.0*l+1.0) * factorial(l-m)/(4.0*M_PI*factorial(l+m))) * 
         legendre_pml(m, l, cos(theta)) * cexp(CMPLX(0, 1) * m * phi);
}

// FIXME: Remove this when everything is implemented.
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

#endif

