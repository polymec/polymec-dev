// Copyright (c) 2012-2016, Jeffrey N. Johnson
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

// A tuple is an N-tuple of elements with some associated comparison 
// functions. Really, the tuple interface manipulates arrays of 
// length N in an object-oriented fashion. One defines the interface for 
// a tuple with elements of a given type using
// DEFINE_TUPLE(tuple_name, element, element_cmp)
//
// Interface for a tuple with datum x defined with 
// DEFINE_TUPLE(x_tuple, x, x_cmp):
//
// x* x_tuple_new(int N) - Creates a new N-tuple with N uninitialized elements.
// void x_tuple_free(x* tuple) - Frees the resources associated with the tuple.
// int x_tuple_length(x* tuple) - Returns the length (N) of the tuple.
// int x_tuple_cmp(x* tuple1, x* tuple2) - Returns -1 if tuple1 < tuple2, 
//                                                  0 if tuple1 == tuple2,
//                                                  1 if tuple1 > tuple2.
// bool x_tuple_equals(x* tuple1, x* tuple2) - Returns true if tuple1 == tuple2, 
//                                                     false otherwise.
// int x_tuple_hash(x* tuple) - Returns a hash index for the given tuple.
// tuple[i] - The ith value in the tuple.

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

// Define some basic tuple types.
DEFINE_TUPLE(int_tuple, int, int_cmp)
DEFINE_TUPLE(real_tuple, real_t, real_cmp)

#endif
