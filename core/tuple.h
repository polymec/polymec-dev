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
static inline element* tuple_name##_new(int N) \
{ \
  ASSERT(N > 0); \
  size_t offset = MAX((sizeof(int) / sizeof(element)), 1); \
  element* tuple = malloc(sizeof(element) * (N + offset)); \
  int* length = (int*)(tuple); \
  *length = N; \
  return tuple + offset; \
} \
static inline void tuple_name##_free(element* tuple) \
{ \
  size_t offset = MAX((sizeof(int) / sizeof(element)), 1); \
  element* t = tuple - offset; \
  free(t); \
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
DEFINE_TUPLE(double_tuple, double, double_cmp)

#endif
