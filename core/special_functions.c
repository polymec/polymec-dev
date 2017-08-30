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

complex_t sph_harmonic_ylm(int l, int m, real_t theta, real_t phi)
{
  return pow(-1, m) * 
         sqrt((2.0*l+1.0) * factorial(l-m)/(4.0*M_PI*factorial(l+m))) * 
         legendre_pml(m, l, cos(theta)) * cexp(CMPLX(0, 1) * m * phi);
}

