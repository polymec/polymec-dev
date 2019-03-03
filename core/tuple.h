// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TUPLE_H
#define POLYMEC_TUPLE_H

#include "core/polymec.h"
#include "core/comparators.h"
#include "core/hash_functions.h"

/// \addtogroup core core
///@{

/// \def DEFINE_TUPLE(tuple_name, element, element_cmp)
/// Defines a tuple for the given element type. A tuple is an N-tuple of elements
/// with some associated comparison functions. Really, the tuple interface
/// manipulates arrays of length N in an object-oriented fashion. The following
/// interface is defined for a set with map_name `x`.
/// * `x* x_tuple_new(int N)` - Creates a new N-tuple with N uninitialized elements.
/// * `void x_tuple_free(x* tuple)` - Frees the resources associated with the tuple.
/// * `int x_tuple_length(x* tuple)` - Returns the length (N) of the tuple.
/// * `x* x_tuple_clone(x* tuple)` - Creates a new copy of the given tuple.
/// * `void x_tuple_copy(x* src, x* dest)` - Copies the elements of src to dest. The tuples must be of equal length.
/// * `int x_tuple_cmp(x* tuple1, x* tuple2)` - Returns -1 if tuple1 < tuple2, 0 if tuple1 == tuple2, 1 if tuple1 > tuple2.
/// * `bool x_tuple_equals(x* tuple1, x* tuple2)` - Returns true if tuple1 == tuple2, false otherwise.
/// * `int x_tuple_hash(x* tuple)` - Returns a hash index for the given tuple.
///
/// Data members for a tuple `tuple`:
/// * `tuple[i]` - The ith value in the tuple.
///
/// \param tuple_name The name of the tuple.
/// \param element The data type stored in the tuple.
/// \param element_cmp A comparator function that accepts two elements and
///                    returns values specified in the `x_tuple_cmp` interface.
#define DEFINE_TUPLE(tuple_name, element, element_cmp) \
\
typedef element tuple_name##_value_t; \
static inline element* tuple_name##_new(int N) \
{ \
  ASSERT(N > 0); \
  size_t offset = MAX((sizeof(int) / sizeof(element)), 1); \
  element* tuple = (element*)polymec_malloc(sizeof(element) * (N + offset)); \
  int* length = (int*)(tuple); \
  *length = N; \
  return tuple + offset; \
} \
\
static inline void tuple_name##_free(element* tuple) \
{ \
  size_t offset = MAX((sizeof(int) / sizeof(element)), 1); \
  element* t = tuple - offset; \
  polymec_free(t); \
} \
\
static inline int tuple_name##_length(element* tuple) \
{ \
  size_t offset = MAX((sizeof(int) / sizeof(element)), 1); \
  int* length = (int*)(tuple - offset); \
  return *length; \
} \
\
static inline void tuple_name##_copy(element* src, element* dest) \
{ \
  size_t offset = MAX((sizeof(int) / sizeof(element)), 1); \
  int src_len = *((int*)(src - offset)); \
  ASSERT(*((int*)(dest - offset)) == src_len); \
  for (int i = 0; i < src_len; ++i) \
    dest[i] = src[i]; \
} \
\
static inline element* tuple_name##_clone(element* tuple) \
{ \
  element* clone = tuple_name##_new(tuple_name##_length(tuple)); \
  tuple_name##_copy(tuple, clone); \
  return clone; \
} \
\
static inline int tuple_name##_cmp(element* tuple1, element* tuple2) \
{ \
  int N1 = tuple_name##_length(tuple1); \
  int N2 = tuple_name##_length(tuple2); \
  int N_min = MIN(N1, N2); \
  ASSERT(N_min > 0); \
  int cmp = element_cmp(tuple1[0], tuple2[0]); \
  int i = 1; \
  while ((cmp == 0) && (i < N_min)) \
  { \
    cmp = element_cmp(tuple1[i], tuple2[i]); \
    ++i; \
  } \
  if (cmp == 0) \
  { \
    cmp = (N_min < N1) ? 1 : (N_min < N2) ? -1 : 0; \
  } \
  return cmp; \
} \
\
static inline bool tuple_name##_equals(element* tuple1, element* tuple2) \
{ \
  return (tuple_name##_cmp(tuple1, tuple2) == 0); \
} \
\
static inline int tuple_name##_hash(element* tuple) \
{ \
  int N = tuple_name##_length(tuple); \
  return djb2_xor_hash((unsigned char*)tuple, N*sizeof(element)); \
} \

///@}

// Define some basic tuple types.
DEFINE_TUPLE(int_tuple, int, int_cmp)
DEFINE_TUPLE(real_tuple, real_t, real_cmp)

#endif
