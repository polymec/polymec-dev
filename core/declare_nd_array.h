// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DECLARE_ND_ARRAY_H
#define POLYMEC_DECLARE_ND_ARRAY_H

#include "core/polymec.h"

// DECLARE_ND_ARRAY
//
// These macros are used to interpret existing 1D "storage" arrays as 
// multi-dimensional arrays using C99's multidimensional array machinery.
// They take the following arguments.
//   type      - The (primitive) data type of the array.
//   array_var - A multidimensional array that refers to data in storage.
//   storage   - A "flat" array with enough storage to store all the data in 
//               the array.
//   dimX      - The extent of the Xth dimension of the array.

// DECLARE_2D_ARRAY(type, array_var, storage, dim1, dim2)
// DECLARE_3D_ARRAY(type, array_var, storage, dim1, dim2, dim3)
// DECLARE_4D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4)
// DECLARE_5D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5)
// DECLARE_6D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6)
// DECLARE_7D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6, dim7)
// 
// We go to 7D just to make Fortran feel less special. :-)
//
#define DECLARE_2D_ARRAY(type, array_var, storage, dim1, dim2) \
type (*array_var)[dim2] = (void*)storage
#define DECLARE_3D_ARRAY(type, array_var, storage, dim1, dim2, dim3) \
type (*array_var)[dim2][dim3] = (void*)storage
#define DECLARE_4D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4) \
type (*array_var)[dim2][dim3][dim4] = (void*)storage
#define DECLARE_5D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5) \
type (*array_var)[dim2][dim3][dim4][dim5] = (void*)storage
#define DECLARE_6D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6) \
type (*array_var)[dim2][dim3][dim4][dim5][dim6] = (void*)storage
#define DECLARE_7D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6, dim7) \
type (*array_var)[dim2][dim3][dim4][dim5][dim6][dim7] = (void*)storage

// ARRAY_INDEX_ND
//
// These macros are used to map multiple array indices onto a flattened index
// into the underlying block of memory occupied by a multidimensional array.
// They prevent the programmer from having to remember how C represents multi-
// dimensional arrays in memory.

// ARRAY_INDEX_2D(dim1, dim2, i, j) -> I
// ARRAY_INDEX_3D(dim1, dim2, dim3, i, j, k) -> I
// ARRAY_INDEX_4D(dim1, dim2, dim3, dim4, i, j, k, l) -> I
// ARRAY_INDEX_5D(dim1, dim2, dim3, dim4, dim5, i, j, k, l, m) -> I
// ARRAY_INDEX_6D(dim1, dim2, dim3, dim4, dim5, dim6, i, j, k, l, m, n) -> I
// ARRAY_INDEX_7D(dim1, dim2, dim3, dim4, dim5, dim6, dim7, i, j, k, l, m, n, p) -> I
//
#define ARRAY_INDEX_2D(dim1, dim2, i, j) (dim2*i + j)
#define ARRAY_INDEX_3D(dim1, dim2, dim3, i, j, k) (dim2*dim3*i + dim2*j + k)
#define ARRAY_INDEX_4D(dim1, dim2, dim3, dim4, i, j, k, l) (dim2*dim3*dim4*i + dim2*dim3*j + dim2*k + l)
#define ARRAY_INDEX_5D(dim1, dim2, dim3, dim4, dim5, i, j, k, l, m) (dim2*dim3*dim4*dim5*i + dim2*dim3*dim4*j + dim2*dim3*k + dim2*l + m)
#define ARRAY_INDEX_6D(dim1, dim2, dim3, dim4, dim5, dim6, i, j, k, l, m, n) (dim2*dim3*dim4*dim5*dim6*i + dim2*dim3*dim4*dim5*j + dim2*dim3*dim4*k + dim2*dim3*l + dim2*m + n)
#define ARRAY_INDEX_7D(dim1, dim2, dim3, dim4, dim5, dim6, dim7, i, j, k, l, m, n, p) (dim2*dim3*dim4*dim5*dim6*dim7*i + dim2*dim3*dim4*dim5*dim6*j + dim2*dim3*dim4*dim5*k + dim2*dim3*dim4*l + dim2*dim3*m + dim2*n + p)

#endif
