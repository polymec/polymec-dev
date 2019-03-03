// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNORDERED_SET_H
#define POLYMEC_UNORDERED_SET_H

#include "core/unordered_map.h"
#include "core/hash_functions.h"
#include "core/comparators.h"
#include "core/tuple.h"

/// \addtogroup core core
///@{

/// \def DEFINE_UNORDERED_SET(set_name, element, hash_func, equals_func)
/// Defines an unordered set for the given element type. The following
/// interface is defined for a set with map_name `x_set`.
/// * `x_set_t* x_set_new()` - Creates a new empty ordered set.
/// * `void x_set_free(x_set_t* set)` - Destroys the set.
/// * `void x_set_clear(x_set_t* set)` - Empties the set.
/// * `bool x_set_contains(x_set_t* set, x datum)` - Returns true if the set contains the datum, false otherwise.
/// * `void x_set_insert(x_set_t* set, x datum)` - Inserts a datum into the set.
/// * `void x_set_insert_with_dtor(x_set_t* set, x datum, x_dtor dtor)` - Inserts a datum into the set with the given destructor.
/// * `void x_set_delete(x_set_t* set, x datum)` - Deletes the datum from the set.
/// * `bool x_set_next(x_set_t* set, int* pos, x* datum)` - Allows traversal of the set.
/// * `void x_set_union(x_set_t* set, x_set_t* other_set, x_set_t* union)` - Differences this set with the other set, storing the result in intersection.
/// * `void x_set_intersection(x_set_t* set, x_set_t* other_set, x_set_t* intersection)` - Intersects this set with the other set, storing the result in intersection.
/// * `void x_set_difference(x_set_t* set, x_set_t* other_set, x_set_t* difference)` - Differences this set with the other set, storing the result in intersection.
/// * `bool x_set_empty(x_set_t* set)` - Returns true if empty, false otherwise.
/// * `set->size` - The size of the set.
///
/// Member data for an unordered set `set`:
/// * `set->size` - The number of elements in the set.
///
/// \param set_name The name of the unordered set.
/// \param element The data type stored by the set.
/// \param hash_func A hash function mapping an element to an integer.
/// \param equals_func A comparator function that accepts two elements and
///                    returns true if these elements are equal, false otherwise.
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
  set_name##_t* set = (set_name##_t*)polymec_malloc(sizeof(set_name##_t)); \
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
  polymec_free(set); \
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
static inline void set_name##_union(set_name##_t* set, set_name##_t* other_set, set_name##_t* union_) \
{ \
  set_name##_clear(union_); \
  int pos = 0; \
  set_name##_element_t e; \
  while (set_name##_next(set, &pos, &e)) \
    set_name##_insert(union_, e); \
  pos = 0; \
  while (set_name##_next(other_set, &pos, &e)) \
    set_name##_insert(union_, e); \
} \
\
static inline void set_name##_intersection(set_name##_t* set, set_name##_t* other_set, set_name##_t* intersection) \
{ \
  set_name##_clear(intersection); \
  int pos = 0; \
  set_name##_element_t e; \
  while (set_name##_next(set, &pos, &e)) \
  { \
    if (set_name##_contains(other_set, e)) \
      set_name##_insert(intersection, e); \
  } \
} \
\
static inline void set_name##_difference(set_name##_t* set, set_name##_t* other_set, set_name##_t* difference) \
{ \
  set_name##_clear(difference); \
  int pos = 0; \
  set_name##_element_t e; \
  while (set_name##_next(set, &pos, &e)) \
  { \
    if (!set_name##_contains(other_set, e)) \
      set_name##_insert(difference, e); \
  } \
} \
\
static inline bool set_name##_empty(set_name##_t* set) \
{ \
  return (set->size == 0); \
} \
\

///@}

// Define some ordered_sets.

/// \class int_unordered_set
/// An unordered set of integers.
DEFINE_UNORDERED_SET(int_unordered_set, int, int_hash, int_equals)

/// \class index_unordered_set
/// An unordered set of 64-bit indices.
DEFINE_UNORDERED_SET(index_unordered_set, index_t, index_hash, index_equals)

/// \class string_unordered_set
/// An unordered set of strings (char*).
DEFINE_UNORDERED_SET(string_unordered_set, char*, string_hash, string_equals)

/// \class int_tuple_unordered_set
/// An unordered set of int tuples (int*).
DEFINE_UNORDERED_SET(int_tuple_unordered_set, int*, int_tuple_hash, int_tuple_equals)

/// \class int_pair_unordered_set
/// An unordered set of int pairs (int[2]).
DEFINE_UNORDERED_SET(int_pair_unordered_set, int*, int_pair_hash, int_pair_equals)

/// \class ptr_unordered_set
/// An unordered set of pointers.
DEFINE_UNORDERED_SET(ptr_unordered_set, void*, ptr_hash, ptr_equals)

#endif
