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

#ifndef POLYMEC_UNORDERED_SET_H
#define POLYMEC_UNORDERED_SET_H

#include "core/unordered_map.h"
#include "core/hash_functions.h"
#include "core/comparators.h"
#include "core/tuple.h"

// An unordered set is a container that stores unique values using a hash table.
// One defines an unordered set using
// DEFINE_UNORDERED_SET(set_name, element, hash_func, equals_func)
//
// Interface for a type x_set_t (with node type x_set_node_t 
// and datum x) defined with 
// DEFINE_UNORDERED_SET(x_set, x, x_hash, x_equals):
// 
// x_set_t* x_set_new() - Creates a new empty ordered set.
// void x_set_free(x_set_t* set) - Destroys the set.
// void x_set_clear(x_set_t* set) - Empties the set.
// bool x_set_contains(x_set_t* set, x datum) - Returns true if the set contains the datum, false otherwise.
// void x_set_insert(x_set_t* set, x datum) - Inserts a datum into the set.
// void x_set_delete(x_set_t* set, x datum) - Deletes the datum from the set.
// void x_set_next(x_set_t* set, int* pos, x* datum) - Allows traversal of the set.

#define DEFINE_UNORDERED_SET(set_name, element, hash_func, equals_func) \
DEFINE_UNORDERED_MAP(set_name##_unordered_map, element, bool, hash_func, equals_func) \
typedef element set_name##_element_t; \
typedef int (*set_name##_hash_func)(element); \
typedef set_name##_unordered_map_k_dtor set_name##_dtor; \
typedef struct \
{ \
  set_name##_unordered_map_t* map; \
  int size; \
} set_name##_t; \
\
static inline set_name##_t* set_name##_new() \
{ \
  set_name##_t* set = malloc(sizeof(set_name##_t)); \
  set->map = set_name##_unordered_map_new(); \
  set->size = 0; \
  return set; \
} \
\
static inline void set_name##_clear(set_name##_t* set) \
{ \
  set_name##_unordered_map_clear(set->map); \
  set->size = 0; \
} \
\
static inline void set_name##_free(set_name##_t* set) \
{ \
  set_name##_unordered_map_free(set->map); \
  free(set); \
} \
\
static inline bool set_name##_contains(set_name##_t* set, set_name##_element_t datum) \
{ \
  return set_name##_unordered_map_contains(set->map, datum); \
} \
\
static inline void set_name##_insert_with_dtor(set_name##_t* set, set_name##_element_t datum, set_name##_dtor dtor) \
{ \
  set_name##_unordered_map_insert_with_k_dtor(set->map, datum, true, dtor); \
  set->size = set->map->size; \
} \
\
static inline void set_name##_insert(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_insert_with_dtor(set, datum, NULL); \
  set->size = set->map->size; \
} \
\
static inline void set_name##_delete(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_unordered_map_delete(set->map, datum); \
  set->size = set->map->size; \
} \
\
static inline bool set_name##_next(set_name##_t* set, int* pos, set_name##_element_t* datum) \
{ \
  bool val; \
  return set_name##_unordered_map_next(set->map, pos, datum, &val); \
} \
\

// Define some ordered_sets.
DEFINE_UNORDERED_SET(int_unordered_set, int, int_hash, int_equals)
DEFINE_UNORDERED_SET(string_unordered_set, char*, string_hash, string_equals)
DEFINE_UNORDERED_SET(int_tuple_unordered_set, int*, int_tuple_hash, int_tuple_equals)


#endif
