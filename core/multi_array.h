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

#ifndef POLYMEC_MULTI_ARRAY_H
#define POLYMEC_MULTI_ARRAY_H

#include "core/polymec.h"

// These functions can be used to allocate storage for multi-dimensional 
// arrays that are allocated with locality in mind. The storage requirements 
// for the arrays are somewhat larger than single-dimensional arrays, as 
// the pointers must be stored as well as the array data, so use caution when 
// using multi-dimensional arrays over large index spaces.
// One defines allocation/deallocation functions for a multi-dimensional 
// array using
// DEFINE_MULTI_ARRAY(prefix, element)
//
// The functions for a multi-dimensional array with prefix x and type x are:
// 
// x** x_array2_new(int dim1, int dim2) - Creates a new 2D array on the heap with the given dimensions.
// void x_array2_free(x** array, int dim1, int dim2) - Frees a 2D array.
// x*** x_array3_new(int dim1, int dim2, int dim3) - Creates a new 3D array on the heap with the given dimensions.
// void x_array3_free(x*** array, int dim1, int dim2, int dim3) - Frees a 3D array.
// x**** x_array4_new(int dim1, int dim2, int dim3, dim4) - Creates a new 4D array on the heap with the given dimensions.
// void x_array4_free(x**** array, int dim1, int dim2, int dim3, int dim4) - Frees a 4D array.
// 
// Note that dimension information for the array is not stored with the array 
// itself, and must be tracked separately.

#define DEFINE_MULTI_ARRAY(prefix, element) \
static inline element** prefix##_array2_new(int dim1, int dim2) \
{ \
  ASSERT(dim1 > 0); \
  ASSERT(dim2 > 0); \
  element** array = polymec_malloc(sizeof(element*) * dim1); \
  array[0] = polymec_malloc(sizeof(element) * dim1 * dim2); \
  for (int i = 1; i < dim1; ++i) \
    array[i] = &array[0][i*dim2]; \
  return array; \
} \
\
static inline void prefix##_array2_free(element** array, int dim1, int dim2) \
{ \
  ASSERT(array != NULL); \
  polymec_free(array[0]); \
  polymec_free(array); \
} \
\
static inline element*** prefix##_array3_new(int dim1, int dim2, int dim3) \
{ \
  ASSERT(dim1 > 0); \
  ASSERT(dim2 > 0); \
  ASSERT(dim3 > 0); \
  element*** array = polymec_malloc(sizeof(element**) * dim1); \
  for (int i = 0; i < dim1; ++i) \
  { \
    array[i] = polymec_malloc(sizeof(element*) * dim2); \
    for (int j = 0; j < dim2; ++j) \
    { \
      if ((i == 0) && (j == 0)) \
        array[0][0] = polymec_malloc(sizeof(element) * dim1 * dim2 * dim3); \
      else \
        array[i][j] = &array[0][0][i*dim2*dim3 + j*dim3]; \
    } \
  } \
  return array; \
} \
\
static inline void prefix##_array3_free(element*** array, int dim1, int dim2, int dim3) \
{ \
  ASSERT(array != NULL); \
  for (int i = 0; i < dim1; ++i) \
  { \
    if (i == 0) \
      polymec_free(array[0][0]); \
    polymec_free(array[i]); \
  } \
  polymec_free(array); \
} \
\
static inline element**** prefix##_array4_new(int dim1, int dim2, int dim3, int dim4) \
{ \
  ASSERT(dim1 > 0); \
  ASSERT(dim2 > 0); \
  ASSERT(dim3 > 0); \
  ASSERT(dim4 > 0); \
  element**** array = polymec_malloc(sizeof(element***) * dim1); \
  for (int i = 0; i < dim1; ++i) \
  { \
    array[i] = polymec_malloc(sizeof(element**) * dim2); \
    for (int j = 0; j < dim2; ++j) \
    { \
      array[i][j] = polymec_malloc(sizeof(element*) * dim3); \
      for (int k = 0; k < dim3; ++k) \
      { \
        if ((i == 0) && (j == 0) && (k == 0)) \
          array[0][0][0] = polymec_malloc(sizeof(element) * dim1 * dim2 * dim3 * dim4); \
        else \
          array[i][j][k] = &array[0][0][0][i*dim2*dim3*dim4 + j*dim3*dim4 + k*dim4]; \
      } \
    } \
  } \
  return array; \
} \
\
static inline void prefix##_array4_free(element**** array, int dim1, int dim2, int dim3, int dim4) \
{ \
  ASSERT(array != NULL); \
  for (int i = 0; i < dim1; ++i) \
  { \
    for (int j = 0; j < dim2; ++j) \
    { \
      if ((i == 0) && (j == 0)) \
        polymec_free(array[0][0][0]); \
      polymec_free(array[i][j]); \
    } \
    polymec_free(array[i]); \
  } \
  polymec_free(array); \
}

// Define some basic array types.
DEFINE_MULTI_ARRAY(int, int)
DEFINE_MULTI_ARRAY(real, real_t)

#endif
