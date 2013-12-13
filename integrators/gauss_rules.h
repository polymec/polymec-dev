// Copyright (c) 2012-2013, Jeffrey N. Johnson
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


#ifndef POLYMEC_GAUSS_RULES_H
#define POLYMEC_GAUSS_RULES_H

// Herein are functions for computing various Gaussian quadratures.

// This function fills the given array (size n+1 each) with the quadrature
// points and weights for the order n Gauss integration rule on the 
// interval (­1, 1).
void get_gauss_points(int n, double* points, double* weights);

// This function fills the given array (size n+1 each) with the quadrature
// points and weights for the order n Gauss-Legendre integration rule on the 
// interval (­1, 1).
void get_gauss_legendre_points(int n, double* points, double* weights);

// This function fills the given array (size n+1 each) with the quadrature
// points and weights for the order n Gauss-Radau integration rule on the 
// interval (­1, 1].
void get_gauss_radau_points(int n, double* points, double* weights);

// This function fills the given array (size n+1 each) with the quadrature
// points and weights for the order n Gauss-Lobatto integration rule on the 
// interval [­1, 1].
void get_gauss_lobatto_points(int n, double* points, double* weights);

#endif
