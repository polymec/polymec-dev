// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
