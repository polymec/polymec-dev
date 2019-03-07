// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DECLARE_ND_ARRAY_H
#define POLYMEC_DECLARE_ND_ARRAY_H

#include "core/polymec.h"

///@{
// These macros are used to interpret existing 1D "storage" arrays as
// multi-dimensional arrays using C99's multidimensional array machinery.
// They take the following arguments.
//   type      - The (primitive) data type of the array.
//   array_var - A multidimensional array that refers to data in storage.
//   storage   - A "flat" array with enough storage to store all the data in
//               the array.
//   dimX      - The extent of the Xth dimension of the array.

/// \def DECLARE_2D_ARRAY
/// Declares a 2D array.
/// \param type The (primitive) data type of the array.
/// \param array_var The declared multidimensional array.
/// \param storage A "flat" array with enough room to fit all the data in array_var.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
#define DECLARE_2D_ARRAY(type, array_var, storage, dim1, dim2) \
type (*array_var)[dim2] = (void*)storage

/// \def DECLARE_3D_ARRAY
/// Declares a 3D array.
/// \param type The (primitive) data type of the array.
/// \param array_var The declared multidimensional array.
/// \param storage A "flat" array with enough room to fit all the data in array_var.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
#define DECLARE_3D_ARRAY(type, array_var, storage, dim1, dim2, dim3) \
type (*array_var)[dim2][dim3] = (void*)storage

/// \def DECLARE_4D_ARRAY
/// Declares a 4D array.
/// \param type The (primitive) data type of the array.
/// \param array_var The declared multidimensional array.
/// \param storage A "flat" array with enough room to fit all the data in array_var.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
#define DECLARE_4D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4) \
type (*array_var)[dim2][dim3][dim4] = (void*)storage

/// \def DECLARE_5D_ARRAY
/// Declares a 5D array.
/// \param type The (primitive) data type of the array.
/// \param array_var The declared multidimensional array.
/// \param storage A "flat" array with enough room to fit all the data in array_var.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param dim5 The fifth dimension of the array's component space.
#define DECLARE_5D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5) \
type (*array_var)[dim2][dim3][dim4][dim5] = (void*)storage

/// \def DECLARE_6D_ARRAY
/// Declares a 6D array.
/// \param type The (primitive) data type of the array.
/// \param array_var The declared multidimensional array.
/// \param storage A "flat" array with enough room to fit all the data in array_var.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param dim5 The fifth dimension of the array's component space.
/// \param dim6 The sixth dimension of the array's component space.
#define DECLARE_6D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6) \
type (*array_var)[dim2][dim3][dim4][dim5][dim6] = (void*)storage

/// \def DECLARE_7D_ARRAY
/// Declares a 7D array.
/// \param type The (primitive) data type of the array.
/// \param array_var The declared multidimensional array.
/// \param storage A "flat" array with enough room to fit all the data in array_var.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param dim5 The fifth dimension of the array's component space.
/// \param dim6 The sixth dimension of the array's component space.
/// \param dim7 The seventh dimension of the array's component space.
#define DECLARE_7D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6, dim7) \
type (*array_var)[dim2][dim3][dim4][dim5][dim6][dim7] = (void*)storage

///@}

///@{
// These macros are used to map multiple array indices onto a flattened index
// into the underlying block of memory occupied by a multidimensional array.
// They prevent the programmer from having to remember how C represents multi-
// dimensional arrays in memory.

/// \def ARRAY_INDEX_2D(dim1, dim2, i, j) -> I
/// Returns the storage offset for a set of components in a multidimensional array.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param i The first component in the array's component space.
/// \param j The second component in the array's component space.
/// \returns an index in a flattened index space that corresponds to the
///          components of a 2D array.
#define ARRAY_INDEX_2D(dim1, dim2, i, j) ((dim2)*(i) + (j))

/// \def ARRAY_INDEX_3D(dim1, dim2, dim3, i, j, k) -> I
/// Returns the storage offset for a set of components in a multidimensional array.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param i The first component in the array's component space.
/// \param j The second component in the array's component space.
/// \param k The third component in the array's component space.
/// \returns an index in a flattened index space that corresponds to the
///          components of a 3D array.
#define ARRAY_INDEX_3D(dim1, dim2, dim3, i, j, k) ((dim2)*(dim3)*(i) + (dim3)*(j) + (k))

/// \def ARRAY_INDEX_4D(dim1, dim2, dim3, dim4, i, j, k, l) -> I
/// Returns the storage offset for a set of components in a multidimensional array.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param i The first component in the array's component space.
/// \param j The second component in the array's component space.
/// \param k The third component in the array's component space.
/// \param l The fourth component in the array's component space.
/// \returns an index in a flattened index space that corresponds to the
///          components of a 4D array.
#define ARRAY_INDEX_4D(dim1, dim2, dim3, dim4, i, j, k, l) ((dim2)*(dim3)*(dim4)*(i) + (dim3)*(dim4)*(j) + (dim4)*(k) + (l))

/// \def ARRAY_INDEX_5D(dim1, dim2, dim3, dim4, dim5, i, j, k, l, m) -> I
/// Returns the storage offset for a set of components in a multidimensional array.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param dim5 The fifth dimension of the array's component space.
/// \param i The first component in the array's component space.
/// \param j The second component in the array's component space.
/// \param k The third component in the array's component space.
/// \param l The fourth component in the array's component space.
/// \param m The fifth component in the array's component space.
/// \returns an index in a flattened index space that corresponds to the
///          components of a 5D array.
#define ARRAY_INDEX_5D(dim1, dim2, dim3, dim4, dim5, i, j, k, l, m) ((dim2)*(dim3)*(dim4)*(dim5)*(i) + (dim3)*(dim4)*(dim5)*(j) + (dim4)*(dim5)*(k) + (dim5)*(l) + (m))

/// \def ARRAY_INDEX_6D(dim1, dim2, dim3, dim4, dim5, dim6, i, j, k, l, m, n) -> I
/// Returns the storage offset for a set of components in a multidimensional array.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param dim5 The fifth dimension of the array's component space.
/// \param dim6 The sixth dimension of the array's component space.
/// \param i The first component in the array's component space.
/// \param j The second component in the array's component space.
/// \param k The third component in the array's component space.
/// \param l The fourth component in the array's component space.
/// \param m The fifth component in the array's component space.
/// \param n The sixth component in the array's component space.
/// \returns an index in a flattened index space that corresponds to the
///          components of a 6D array.
#define ARRAY_INDEX_6D(dim1, dim2, dim3, dim4, dim5, dim6, i, j, k, l, m, n) ((dim2)*(dim3)*(dim4)*(dim5)*(dim6)*(i) + (dim3)*(dim4)*(dim5)*(dim6)*(j) + (dim4)*(dim5)*(dim6)*(k) + (dim5)*(dim6)*(l) + (dim6)*(m) + (n))

/// \def ARRAY_INDEX_7D(dim1, dim2, dim3, dim4, dim5, dim6, dim7, i, j, k, l, m, n, p) -> I
/// Returns the storage offset for a set of components in a multidimensional array.
/// \param dim1 The first dimension of the array's component space.
/// \param dim2 The second dimension of the array's component space.
/// \param dim3 The third dimension of the array's component space.
/// \param dim4 The fourth dimension of the array's component space.
/// \param dim5 The fifth dimension of the array's component space.
/// \param dim6 The sixth dimension of the array's component space.
/// \param dim7 The seventh dimension of the array's component space.
/// \param i The first component in the array's component space.
/// \param j The second component in the array's component space.
/// \param k The third component in the array's component space.
/// \param l The fourth component in the array's component space.
/// \param m The fifth component in the array's component space.
/// \param n The sixth component in the array's component space.
/// \param p The seventh component in the array's component space.
/// \returns an index in a flattened index space that corresponds to the
///          components of a 7D array.
#define ARRAY_INDEX_7D(dim1, dim2, dim3, dim4, dim5, dim6, dim7, i, j, k, l, m, n, p) ((dim2)*(dim3)*(dim4)*(dim5)*(dim6)*(dim7)*(i) + (dim3)*(dim4)*(dim5)*(dim6)*(dim7)*(j) + (dim4)*(dim5)*(dim6)*(dim7)*(k) + (dim5)*(dim6)*(dim7)*(l) + (dim6)*(dim7)*(m) + (dim7)*(n) + (p))

///@}

#endif
