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
  element* tuple = polymec_malloc(sizeof(element) * (N + offset)); \
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
  ASSERT(N1 > 0); \
  ASSERT(tuple_name##_length(tuple2) > 0); \
  ASSERT(N1 == tuple_name##_length(tuple2)); \
  int cmp = element_cmp(tuple1[0], tuple2[0]); \
  int i = 1; \
  while ((cmp == 0) && (i < N1)) \
  { \
    cmp = element_cmp(tuple1[i], tuple2[i]); \
    ++i; \
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
