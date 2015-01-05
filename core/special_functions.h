// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SPECIAL_FUNCTIONS_H
#define POLYMEC_SPECIAL_FUNCTIONS_H

// These functions implement special functions not available in the 
// standard C math library. They are adapted from the Fortran 77 subroutines 
// that accompany "Computation of Special Functions" by Shanjie Zhang and 
// Jianming Jin (Copyright (c) 1996, John Wiley and Sons, Inc).
// This copyright/attribution must accompany any reproduction of these 
// functions.

// Returns the gamma function with the given real argument. The argument 
// cannot exceed 171.6, nor can it be a negative integer.
double gamma(double x);

// Returns the value of J0, the Bessel function of the first kind of order 0, 
// at the given real value x.
double bessel_j0(double x);

// Returns the value of J1, the Bessel function of the first kind of order 1, 
// at the given real value x.
double bessel_j1(double x);

// Returns the value of Jn, the Bessel function of the first kind of order n, 
// at the given real value x.
double bessel_jn(int n, double x);

// Returns the value of the first derivative of J0, the Bessel function of 
// the first kind of order 0, at the given real value x.
double bessel_dj0dx(double x);

// Returns the value of the first derivative of J1, the Bessel function of 
// the first kind of order 1, at the given real value x.
double bessel_dj1dx(double x);

// Returns the value of the first derivative of Jn, the Bessel function of 
// the first kind of order n, at the given real value x.
double bessel_djndx(int n, double x);

// Populates the given array with the specified number of positive roots of 
// Jn, the Bessel function of the first kind of order n. The roots array must 
// be able to store the requested roots. 
void bessel_find_jn_roots(int n, int num_roots, double* roots);

// Returns the value of Y0, the Bessel function of the second kind of order 0, 
// at the given real value x.
double bessel_y0(double x);

// Returns the value of Y1, the Bessel function of the second kind of order 1, 
// at the given real value x.
double bessel_y1(double x);

// Returns the value of Yn, the Bessel function of the second kind of order n, 
// at the given real value x.
double bessel_yn(int n, double x);

// Returns the value of the first derivative of Y0, the Bessel function of 
// the second kind of order 0, at the given real value x.
double bessel_dy0dx(double x);

// Returns the value of the first derivative of Y1, the Bessel function of 
// the second kind of order 1, at the given real value x.
double bessel_dy1dx(double x);

// Returns the value of the first derivative of Yn, the Bessel function of 
// the second kind of order n, at the given real value x.
double bessel_dyndx(int n, double x);

// Populates the given array with the specified number of positive roots of 
// Yn, the Bessel function of the first second of order n. The roots array 
// must be able to store the requested roots. 
void bessel_find_yn_roots(int n, int num_roots, double* roots);

// Returns the value of I0, the modified Bessel function of the first kind of 
// order zero, at the given real value x.
double bessel_i0(double x);

// Returns the value of I1, the modified Bessel function of the first kind of 
// order one, at the given real value x.
double bessel_i1(double x);

// Returns the value of In, the modified Bessel function of the first kind of 
// order n, at the given real value x.
double bessel_in(int n, double x);

// Returns the value of the first derivative of I0, the modified Bessel function of 
// the first kind of order 0, at the given real value x.
double bessel_di0dx(double x);

// Returns the value of the first derivative of I1, the modified Bessel function of 
// the first kind of order 1, at the given real value x.
double bessel_di1dx(double x);

// Returns the value of the first derivative of In, the modified Bessel function of 
// the first kind of order n, at the given real value x.
double bessel_dindx(int n, double x);

// Returns the value of K0, the modified Bessel function of the second kind of 
// order zero, at the given real value x.
double bessel_k0(double x);

// Returns the value of K1, the modified Bessel function of the second kind of 
// order one, at the given real value x.
double bessel_k1(double x);

// Returns the value of Kn, the modified Bessel function of the second kind of 
// order n, at the given real value x.
double bessel_kn(int n, double x);

// Returns the value of the first derivative of K0, the modified Bessel function of 
// the second kind of order 0, at the given real value x.
double bessel_dk0dx(double x);

// Returns the value of the first derivative of K1, the modified Bessel function of 
// the second kind of order 1, at the given real value x.
double bessel_dk1dx(double x);

// Returns the value of the first derivative of Kn, the modified Bessel function of 
// the second kind of order n, at the given real value x.
double bessel_dkndx(int n, double x);

// C++ doesn't have complex datatypes, so we don't expose them to C++ compilers.
#ifndef __cplusplus
#include <complex.h>

// Returns the value of J0, the Bessel function of the first kind of order 
// zero, at the given complex value z.
double complex bessel_cj0(double complex z);

// Returns the value of J1, the Bessel function of the first kind of order 
// one, at the given complex value z.
double complex bessel_cj1(double complex z);

// Returns the value of Jn, the Bessel function of the first kind of order n, 
// at the given complex value z.
double complex bessel_cjn(int n, double complex z);

// Returns the value of Y0, the Bessel function of the second kind of order 
// zero, at the given complex value z.
double complex bessel_cy0(double complex z);

// Returns the value of Y1, the Bessel function of the second kind of order
// one, at the given complex value z.
double complex bessel_cy1(double complex z);

// Returns the value of Yn, the Bessel function of the second kind of order n,
// at the given complex value z.
double complex bessel_cyn(int n, double complex z);

// Returns the value of I0, the modified Bessel function of the first kind of 
// order zero, at the given complex value z.
double complex bessel_ci0(double complex z);

// Returns the value of I1, the modified Bessel function of the first kind of 
// order one, at the given complex value z.
double complex bessel_ci1(double complex z);

// Returns the value of In, the modified Bessel function of the first kind of 
// order n, at the given complex value z.
double complex bessel_cin(int n, double complex z);

// Returns the value of K0, the modified Bessel function of the second kind of 
// order zero, at the given complex value z.
double complex bessel_ck0(double complex z);

// Returns the value of K1, the modified Bessel function of the second kind of 
// order one, at the given complex value z.
double complex bessel_ck1(double complex z);

// Returns the value of Kn, the modified Bessel function of the second kind of 
// order n, at the given complex value z.
double complex bessel_ckn(int n, double complex z);

#endif

#endif

