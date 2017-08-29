// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SPECIAL_FUNCTIONS_H
#define POLYMEC_SPECIAL_FUNCTIONS_H

#include "core/polymec.h"

// These functions implement special functions not available in the 
// standard C math library. Several of them are adapted from the Fortran 77 
// subroutines that accompany "Computation of Special Functions" by Shanjie 
// Zhang and Jianming Jin (Copyright (c) 1996, John Wiley and Sons, Inc).
// This copyright/attribution must accompany any reproduction of these 
// functions.

// Returns the value of J0, the Bessel function of the first kind of order 0, 
// at the given real value x.
real_t bessel_j0(real_t x);

// Returns the value of J1, the Bessel function of the first kind of order 1, 
// at the given real value x.
real_t bessel_j1(real_t x);

// Returns the value of Jn, the Bessel function of the first kind of order n, 
// at the given real value x.
real_t bessel_jn(int n, real_t x);

// Returns the value of the first derivative of J0, the Bessel function of 
// the first kind of order 0, at the given real value x.
real_t bessel_dj0dx(real_t x);

// Returns the value of the first derivative of J1, the Bessel function of 
// the first kind of order 1, at the given real value x.
real_t bessel_dj1dx(real_t x);

// Returns the value of the first derivative of Jn, the Bessel function of 
// the first kind of order n, at the given real value x.
real_t bessel_djndx(int n, real_t x);

// Populates the given array with the specified number of positive roots of 
// Jn, the Bessel function of the first kind of order n. The roots array must 
// be able to store the requested roots. 
void bessel_find_jn_roots(int n, 
                          int num_roots, 
                          real_t* roots);

// Returns the value of Y0, the Bessel function of the second kind of order 0, 
// at the given real value x.
real_t bessel_y0(real_t x);

// Returns the value of Y1, the Bessel function of the second kind of order 1, 
// at the given real value x.
real_t bessel_y1(real_t x);

// Returns the value of Yn, the Bessel function of the second kind of order n, 
// at the given real value x.
real_t bessel_yn(int n, real_t x);

// Returns the value of the first derivative of Y0, the Bessel function of 
// the second kind of order 0, at the given real value x.
real_t bessel_dy0dx(real_t x);

// Returns the value of the first derivative of Y1, the Bessel function of 
// the second kind of order 1, at the given real value x.
real_t bessel_dy1dx(real_t x);

// Returns the value of the first derivative of Yn, the Bessel function of 
// the second kind of order n, at the given real value x.
real_t bessel_dyndx(int n, real_t x);

// Populates the given array with the specified number of positive roots of 
// Yn, the Bessel function of the first second of order n. The roots array 
// must be able to store the requested roots. 
void bessel_find_yn_roots(int n, 
                          int num_roots, 
                          real_t* roots);

// Returns the value of I0, the modified Bessel function of the first kind of 
// order zero, at the given real value x.
real_t bessel_i0(real_t x);

// Returns the value of I1, the modified Bessel function of the first kind of 
// order one, at the given real value x.
real_t bessel_i1(real_t x);

// Returns the value of In, the modified Bessel function of the first kind of 
// order n, at the given real value x.
real_t bessel_in(int n, real_t x);

// Returns the value of the first derivative of I0, the modified Bessel function of 
// the first kind of order 0, at the given real value x.
real_t bessel_di0dx(real_t x);

// Returns the value of the first derivative of I1, the modified Bessel function of 
// the first kind of order 1, at the given real value x.
real_t bessel_di1dx(real_t x);

// Returns the value of the first derivative of In, the modified Bessel function of 
// the first kind of order n, at the given real value x.
real_t bessel_dindx(int n, real_t x);

// Returns the value of K0, the modified Bessel function of the second kind of 
// order zero, at the given real value x.
real_t bessel_k0(real_t x);

// Returns the value of K1, the modified Bessel function of the second kind of 
// order one, at the given real value x.
real_t bessel_k1(real_t x);

// Returns the value of Kn, the modified Bessel function of the second kind of 
// order n, at the given real value x.
real_t bessel_kn(int n, real_t x);

// Returns the value of the first derivative of K0, the modified Bessel function of 
// the second kind of order 0, at the given real value x.
real_t bessel_dk0dx(real_t x);

// Returns the value of the first derivative of K1, the modified Bessel function of 
// the second kind of order 1, at the given real value x.
real_t bessel_dk1dx(real_t x);

// Returns the value of the first derivative of Kn, the modified Bessel function of 
// the second kind of order n, at the given real value x.
real_t bessel_dkndx(int n, real_t x);

// Returns the value of Tn, the Chebyshev polynomial of order n, at the given 
// value of x.
real_t chebyshev_tn(int n, real_t x);

// Returns the value of the first derivative of Tn, the Chebyshev polynomial
// of order n, at the given value of x.
real_t chebyshev_dtndx(int n, real_t x);

// Returns the value of Un, the Chebyshev polynomial of order n, at the given 
// value of x.
real_t chebyshev_un(int n, real_t x);

// Returns the value of the first derivative of Un, the Chebyshev polynomial 
// of order n, at the given value of x.
real_t chebyshev_dundx(int n, real_t x);

// Returns the value of Ln, the Laguerre polynomial of order n, at the given 
// value of x.
real_t laguerre_ln(int n, real_t x);

// Returns the value of the first derivative of Ln, the Laguerre polynomial 
// of order n, at the given value of x.
real_t laguerre_dlndx(int n, real_t x);

// Returns the value of Hn, the Hermite polynomial of order n, at the given 
// value of x.
real_t hermite_hn(int n, real_t x);

// Returns the value of dHndx, the first derivative of the Hermite polynomial 
// of order n, at the given value of x.
real_t hermite_dhndx(int n, real_t x);

// Returns the value of Pn, the Legendre polynomial of degree n, at the given 
// real value of x.
real_t legendre_pn(int n, real_t x);

// Returns the value of Pml, the associated Legendre polynomial of order m 
// and degree l, at the given real value x.
real_t legendre_pml(int m, int l, real_t x);

// C++ doesn't have complex datatypes, so we don't expose complex functions 
// to those compilers.
#ifndef __cplusplus

// Returns the value of H0, the Bessel function of the third kind (Hankel
// function) of order zero, at the given real value x.
complex_t bessel_h0(real_t x);

// Returns the value of H1, the Bessel function of the third kind (Hankel 
// function) of order one, at the given real value x.
complex_t bessel_h1(real_t x);

// Returns the value of Hn, the Bessel function of the third kind (Hankel
// function) of order n, at the given real value x.
complex_t bessel_hn(int n, real_t x);

// Returns the value of the first derivative of H0, the Bessel function of 
// the third kind (Hankel function) of order 0, at the given real value x.
complex_t bessel_dh0dx(real_t x);

// Returns the value of the first derivative of H1, the Bessel function of 
// the third kind (Hankel function) of order 1, at the given real value x.
complex_t bessel_dh1dx(real_t x);

// Returns the value of the first derivative of Hn, the Bessel function of 
// the third kind (Hankel function) of order n, at the given real value x.
complex_t bessel_dhndx(int n, real_t x);

// Returns the value of J0, the Bessel function of the first kind of order 
// zero, at the given complex value z.
complex_t bessel_cj0(complex_t z);

// Returns the value of J1, the Bessel function of the first kind of order 
// one, at the given complex value z.
complex_t bessel_cj1(complex_t z);

// Returns the value of Jn, the Bessel function of the first kind of order n, 
// at the given complex value z.
complex_t bessel_cjn(int n, complex_t z);

// Returns the value of Jv, the Bessel function of the first kind of 
// complex order nu (written v), at the given complex value z.
complex_t bessel_jv(complex_t v, complex_t z);

// Returns the value of the first derivative of J0, the Bessel function of 
// the first kind of order 0, at the given complex value z.
complex_t bessel_dj0dz(complex_t z);

// Returns the value of the first derivative of J1, the Bessel function of 
// the first kind of order 1, at the given complex value z.
complex_t bessel_dj1dz(complex_t z);

// Returns the value of the first derivative of Jn, the Bessel function of 
// the first kind of integer order n, at the given complex value z.
complex_t bessel_djndz(int n, complex_t z);

// Returns the value of the first derivative of Jv, the Bessel function of 
// the first kind of complex order nu (written v), at the given complex 
// value z.
complex_t bessel_djvdz(complex_t v, complex_t z);

// Returns the value of Y0, the Bessel function of the second kind of order 
// zero, at the given complex value z.
complex_t bessel_cy0(complex_t z);

// Returns the value of Y1, the Bessel function of the second kind of order
// one, at the given complex value z.
complex_t bessel_cy1(complex_t z);

// Returns the value of Yn, the Bessel function of the second kind of order n,
// at the given complex value z.
complex_t bessel_cyn(int n, complex_t z);

// Returns the value of Yv, the Bessel function of the second kind of 
// complex order nu (written v), at the given complex value z.
complex_t bessel_yv(complex_t v, complex_t z);

// Returns the value of the first derivative of Y0, the Bessel function of 
// the second kind of order 0, at the given complex value z.
complex_t bessel_dy0dz(complex_t z);

// Returns the value of the first derivative of Y1, the Bessel function of 
// the second kind of order 1, at the given complex value z.
complex_t bessel_dy1dz(complex_t z);

// Returns the value of the first derivative of Yn, the Bessel function of 
// the second kind of integer order n, at the given complex value z.
complex_t bessel_dyndz(int n, complex_t z);

// Returns the value of the first derivative of Yv, the Bessel function of 
// the second kind of complex order nu (written v), at the given complex 
// value z.
complex_t bessel_dyvdz(complex_t v, complex_t z);

// Returns the value of I0, the modified Bessel function of the first kind of 
// order zero, at the given complex value z.
complex_t bessel_ci0(complex_t z);

// Returns the value of I1, the modified Bessel function of the first kind of 
// order one, at the given complex value z.
complex_t bessel_ci1(complex_t z);

// Returns the value of In, the modified Bessel function of the first kind of 
// order n, at the given complex value z.
complex_t bessel_cin(int n, complex_t z);

// Returns the value of Iv, the modified Bessel function of the first kind of 
// complex order nu (written v), at the given complex value z.
complex_t bessel_iv(complex_t v, complex_t z);

// Returns the value of the first derivative of I0, the modified Bessel 
// function of the first kind of order 0, at the given complex value z.
complex_t bessel_di0dz(complex_t z);

// Returns the value of the first derivative of I1, the modified Bessel 
// function of the first kind of order 1, at the given complex value z.
complex_t bessel_di1dz(complex_t z);

// Returns the value of the first derivative of In, the modified Bessel 
// function of the first kind of integer order n, at the given complex value z.
complex_t bessel_dindz(int n, complex_t z);

// Returns the value of the first derivative of Iv, the modified Bessel 
// function of the first kind of complex order nu (written v), at the given 
// complex value z.
complex_t bessel_divdz(complex_t v, complex_t z);

// Returns the value of K0, the modified Bessel function of the second kind of 
// order zero, at the given complex value z.
complex_t bessel_ck0(complex_t z);

// Returns the value of K1, the modified Bessel function of the second kind of 
// order one, at the given complex value z.
complex_t bessel_ck1(complex_t z);

// Returns the value of Kn, the modified Bessel function of the second kind of 
// order n, at the given complex value z.
complex_t bessel_ckn(int n, complex_t z);

// Returns the value of Kv, the modified Bessel function of the second kind of 
// complex order v (written v), at the given complex value z.
complex_t bessel_kv(int n, complex_t z);

// Returns the value of the first derivative of K0, the modified Bessel 
// function of the second kind of order 0, at the given complex value z.
complex_t bessel_dk0dz(complex_t z);

// Returns the value of the first derivative of K1, the modified Bessel 
// function of the second kind of order 1, at the given complex value z.
complex_t bessel_dk1dz(complex_t z);

// Returns the value of the first derivative of Kn, the modified Bessel 
// function of the second kind of integer order n, at the given complex value z.
complex_t bessel_dkndz(int n, complex_t z);

// Returns the value of the first derivative of Kv, the modified Bessel 
// function of the second kind of complex order nu (written v), at the given 
// complex value z.
complex_t bessel_dkvdz(complex_t v, complex_t z);

// Returns the (complex) value of Ylm, the normalized spherical harmonic, for 
// which -l <= m <= +l, for the real colatitude theta and azimuth phi.
complex_t sph_harmonic_ylm(int l, int m, real_t theta, real_t phi);

#endif

#endif

