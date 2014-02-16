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

#ifndef POLYMEC_COUPLED_LEAST_SQUARES_H
#define POLYMEC_COUPLED_LEAST_SQUARES_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/polynomial.h"

// This file contains logic that implements the Coupled Least Squares 
// algorithm for the construction of compact-stencil high-order polynomial
// fits, as discussed by Haider (2011).
//
// The following concepts are used in this algorithm:
// - A neighborhood of cells W(i) is a set of cells on which a least squares 
//   fit is performed for a quantity u in the vicinity of a cell i. The 
//   neighborhood V(i) refers to the set of cells adjacent to and including i.
// - The moments z(k), associated with a neighborhood W, are volume integrals 
//   involving the kth-order tensor products of differences between the 
//   integration coordinates and the cell center coordinates. Within a given 
//   cell, z(k) is a symmetric tensor of degree k.
// - The linear map w(m|k), given a vector of N cell-averaged values of u on 
//   a neighborhood W, produces the components of the mth derivative of the 
//   degree-k polynomial representation of u within that neighborhood. These 
//   derivatives form a symmetric tensor of rank m, and the tensor contains 
//   6**m components in 3D. Thus, the map w(m|k) can be represented by a 
//   6**m x N matrix. The product of this matrix with the N cell-averaged 
//   values of u in the neighborhood W(i) about a cell i is a vector containing 
//   the 6**m spatial derivatives of u's polynomial representation in i.
//
// Given the linear map w(k|k) for the neighborhood W(i) and the moments 
// z(k+1) for each cell in that neighborhood, this function calculates the 
// components of the linear map w(k+1|k+1). Arguments:
// k - the degree of the derivatives used to construct the k+1 derivatives.
// N - the number of cells in the neighborhood W(i).
// wkk - the components of the linear map w(k|k) in column-major order.
//       There are 6**m * N components in this array.
// zk1_moments - An array containing the N moment tensors z(k+1) for the cells 
//               in the neighborhood W(i), stored in cell-major order. Each 
//               tensor in the array is stored in column-major order.
//               There are N * 3**(k+1) components in this array.
// wk1k1 - an array that will store the components of the linear map w(k+1|k+1)
//         in column-major order.
void reconstruct_cls_derivatives(int k, 
                                 int N,
                                 real_t* wkk, 
                                 real_t* zk1_moments,
                                 real_t* wk1k1);

#endif
