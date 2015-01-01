// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_SPECIAL_FUNCTIONS_H
#define POLYMEC_SPECIAL_FUNCTIONS_H

#include <math.h>

// These functions augment the special functions available in the standard C
// math library. They are adapted from the Fortran 77 subroutines that 
// accompany "Computation of Special Functions" by Shanjie Zhang and 
// Jianming Jin (Copyright (c) 1996, John Wiley and Sons, Inc).

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
// order zero at the given real value x, K0(x).
double k0(double x);

// Returns the value of the modified Bessel function of the second kind, of 
// order one at the given real value x, K1(x).
double k1(double x);

// Returns the value of the modified Bessel function of the second kind, of 
// order n at the given real value x, Kn(x).
double kn(int n, double x);

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
// kind, of order zero at the given complex value x, K0(x).
double complex ck0(double complex x);

// Returns the value of the complex modified Bessel function of the second 
// kind, of order one at the given complex value x, K1(x).
double complex ck1(double complex x);

// Returns the value of the complex modified Bessel function of the second 
// kind, of order n at the given complex value x, Kn(x).
double complex ckn(int n, double complex x);

#endif

#endif

