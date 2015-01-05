// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SPECIAL_FUNCTIONS_H
#define POLYMEC_SPECIAL_FUNCTIONS_H

// These functions augment the special functions available in the standard C
// math library. They are adapted from the Fortran 77 subroutines that 
// accompany "Computation of Special Functions" by Shanjie Zhang and 
// Jianming Jin (Copyright (c) 1996, John Wiley and Sons, Inc).
// This copyright/attribution must accompany any reproduction of these 
// functions.

// Returns the gamma function with the given real argument. The argument 
// cannot exceed 171.6, nor can it be a negative integer.
double gamma(double x);

// Returns the value of the modified Bessel function of the first kind, of 
// order zero at the given real value x: I0(x).
double i0(double x);

// Returns the value of the modified Bessel function of the first kind, of 
// order one at the given real value x: I1(x).
double i1(double x);

// Returns the value of the modified Bessel function of the first kind, of 
// order n at the given real value x: In(x).
double in(int n, double x);

// Returns the value of the modified Bessel function of the second kind, of 
// order zero at the given real value x: K0(x).
double k0(double x);

// Returns the value of the modified Bessel function of the second kind, of 
// order one at the given real value x: K1(x).
double k1(double x);

// Returns the value of the modified Bessel function of the second kind, of 
// order n at the given real value x: Kn(x).
double kn(int n, double x);

// Returns the value of the first derivative of j0, the Bessel function of 
// the first kind of order 0, at the given real value x.
double dj0dx(double x);

// Returns the value of the first derivative of j1, the Bessel function of 
// the first kind of order 1, at the given real value x.
double dj1dx(double x);

// Returns the value of the first derivative of y0, the Bessel function of 
// the second kind of order 0, at the given real value x.
double dy0dx(double x);

// Returns the value of the first derivative of y1, the Bessel function of 
// the second kind of order 1, at the given real value x.
double dy1dx(double x);

// Returns the value of the first derivative of jn, the Bessel function of 
// the first kind of order n, at the given real value x.
double djndx(int n, double x);

// Returns the value of the first derivative of yn, the Bessel function of 
// the second kind of order n, at the given real value x.
double dyndx(int n, double x);

// Populates the given array with the specified number of positive roots of 
// jn, the Bessel function of the first kind of order n. The roots array must 
// be able to store the requested roots. 
void find_jn_roots(int n, int num_roots, double* roots);

// Populates the given array with the specified number of positive roots of 
// yn, the Bessel function of the first second of order n. The roots array 
// must be able to store the requested roots. 
void find_yn_roots(int n, int num_roots, double* roots);

// C++ doesn't have complex datatypes, so we don't expose them to C++ compilers.
#ifndef __cplusplus
#include <complex.h>

// Returns the value of the complex Bessel function of the first kind, 
// of order zero at the given complex value z: J0(z).
double complex cj0(double complex z);

// Returns the value of the complex Bessel function of the first kind, 
// of order one at the given complex value z: J1(z).
double complex cj1(double complex z);

// Returns the value of the complex Bessel function of the first kind, 
// of order n at the given complex value z: Jn(z).
double complex cjn(int n, double complex z);

// Returns the value of the complex Bessel function of the second kind, 
// of order zero at the given complex value z: Y0(z).
double complex cy0(double complex z);

// Returns the value of the complex Bessel function of the second kind, 
// of order one at the given complex value z: Y1(z).
double complex cy1(double complex z);

// Returns the value of the complex Bessel function of the second kind, 
// of order n at the given complex value z: Yn(z).
double complex cyn(int n, double complex z);

// Returns the value of the complex modified Bessel function of the first 
// kind, of order zero at the given complex value z: I0(z).
double complex ci0(double complex z);

// Returns the value of the complex modified Bessel function of the first 
// kind, of order one at the given complex value z: I1(z).
double complex ci1(double complex z);

// Returns the value of the complex modified Bessel function of the first 
// kind, of order n at the given complex value z: In(z).
double complex cin(int n, double complex z);

// Returns the value of the complex modified Bessel function of the second 
// kind, of order zero at the given complex value z: K0(z).
double complex ck0(double complex z);

// Returns the value of the complex modified Bessel function of the second 
// kind, of order one at the given complex value z: K1(z).
double complex ck1(double complex z);

// Returns the value of the complex modified Bessel function of the second 
// kind, of order n at the given complex value z: Kn(z).
double complex ckn(int n, double complex z);

#endif

#endif

