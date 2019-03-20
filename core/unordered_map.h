// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNORDERED_MAP_H
#define POLYMEC_UNORDERED_MAP_H

#include "core/polymec.h"
#include "core/hash_functions.h"
#include "core/comparators.h"
#include "core/tuple.h"

/// \addtogroup core core
///@{

/// \def DEFINE_UNORDERED_MAP(map_name, key_type, value_type, hash_func, equals_func)
/// Defines an unordered map for the given key and value types. The following
/// interface is defined for a map with map_name `x_map`.
/// * `x_map_t* x_map_new(void)` - Creates a new empty unordered map.
/// * `x_map_t* x_map_new_with_capacity(int N)` - Creates a new empty hash map with initial capacity N.
/// * `x_map_t* x_map_clone(x_map_t* m, x_map_key_t (*key_clone_func)(x_map_key_t), x_map_value_t (*val_clone_func)(x_map_value_t))` - Creates a deep copy of m, using the given clone functions to clone the keys and values.
/// * `void x_map_free(x_map_t* map)` - Destroys the map.
/// * `void x_map_clear(x_map_t* map)` - Empties the map.
/// * `x_map_value_t* x_map_get(x_map_t* map, x_map_key_t key)` - Returns the value for the key, or NULL.
/// * `bool x_map_contains(x_map_t* map, x_map_key_t key)` - Returns true if the map contains the key, false if not.
/// * `void x_map_insert(x_map_t* map, x_map_key_t key, x_map_value_t value)` - Sets the value for the given key.
/// * `void x_map_insert_with_kv_dtor(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_kv_dtor)` - Sets the value for the given key and associates a destructor for the key/value pair.
/// * `void x_map_insert_with_kv_dtors(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_k_dtor, x_map_v_dtor)` - Sets the value for the given key and associates separate destructors for the key and value.
/// * `void x_map_insert_with_k_dtor(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_k_dtor)` - Sets the value for the given key and associates a destructor for the key.
/// * `void x_map_insert_with_v_dtor(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_v_dtor)` - Sets the value for the given key and associates a destructor for the value.
/// * `key_type x_map_change_key(x_map_t* map, x_map_key_t old_key, x_map_key_t new_key)` - Renames old_key to new_key, overwriting new_key if it exists. Returns old key.
/// * `void x_map_swap(x_map_t* map, x_map_key_t key1, x_map_key_t key2)` - Swaps the values for key1 and key2, including destructors.
/// * `void x_map_delete(x_map_t* map, x_map_key_t key)` - Deletes the value for the given key.
/// * `bool x_map_next(x_map_t* map, int* pos, x_map_key_t* key, x_map_value_t* value)` - Allows the traversal of the maps keys and values.
/// * `bool x_map_empty(x_map_t* map)` - Returns true if empty, false otherwise.
///
/// Member data for an unordered map `map`:
/// * `map->size` - The number of key-value pairs in the map.
///
/// \param map_name The name of the unordered map.
/// \param key_type The data type used as a key in the map.
/// \param value_type The data type used as a value in the map.
/// \param hash_func A hash function mapping a key to an integer.
/// \param equals_func A comparator function that accepts two key arguments and
///                    returns true if these arguments are equal, false otherwise.

#define DEFINE_UNORDERED_MAP(map_name, key_type, value_type, hash_func, equals_func) \
typedef key_type map_name##_key_t; \
typedef value_type map_name##_value_t; \
typedef int (*map_name##_hash_func)(map_name##_key_t); \
typedef bool (*map_name##_equals_func)(map_name##_key_t, map_name##_key_t); \
typedef void (*map_name##_visitor)(map_name##_key_t, map_name##_value_t, void*); \
typedef struct map_name##_entry_t map_name##_entry_t; \
typedef void (*map_name##_kv_dtor)(key_type, value_type); \
typedef void (*map_name##_k_dtor)(key_type); \
typedef void (*map_name##_v_dtor)(value_type); \
struct map_name##_entry_t \
{ \
  key_type key; \
  int hash; \
  value_type value; \
  map_name##_kv_dtor kv_dtor; \
  map_name##_k_dtor k_dtor; \
  map_name##_v_dtor v_dtor; \
  map_name##_entry_t* next; \
}; \
\
typedef struct \
{ \
  map_name##_entry_t** buckets; \
  int bucket_count; \
  map_name##_hash_func hash; \
  map_name##_equals_func equals; \
  int size; \
  int max_depth; \
} map_name##_t; \
\
static inline map_name##_t* map_name##_new_with_capacity(int N) \
{ \
  map_name##_t* map = (map_name##_t*)polymec_malloc(sizeof(map_name##_t)); \
  int minimum_bucket_count = N * 4 / 3; \
  map->bucket_count = 1; \
  while (map->bucket_count <= minimum_bucket_count) \
    map->bucket_count <<= 1; \
  map->buckets = (map_name##_entry_t**)polymec_calloc(map->bucket_count, sizeof(map_name##_entry_t*)); \
  ASSERT(map->buckets != NULL); \
  map->size = 0; \
  map->hash = hash_func; \
  map->equals = equals_func; \
  map->max_depth = 1; \
  return map; \
} \
\
static inline map_name##_t* map_name##_new(void) \
{ \
  return map_name##_new_with_capacity(32); \
} \
\
static inline void map_name##_clear(map_name##_t* map) \
{ \
  for (int i = 0; i < map->bucket_count; ++i) \
  { \
    map_name##_entry_t* entry = map->buckets[i]; \
    while (entry != NULL) \
    { \
      map_name##_entry_t* next = entry->next; \
      if (entry->kv_dtor != NULL) \
        (*entry->kv_dtor)(entry->key, entry->value); \
      else \
      { \
        if (entry->k_dtor != NULL) \
          (*entry->k_dtor)(entry->key); \
        if (entry->v_dtor != NULL) \
          (*entry->v_dtor)(entry->value); \
      } \
      polymec_free(entry); \
      entry = next; \
    } \
    map->buckets[i] = NULL; \
  } \
  map->size = 0; \
  map->max_depth = 1; \
} \
\
static inline void map_name##_free(map_name##_t* map) \
{ \
  map_name##_clear(map); \
  polymec_free(map->buckets); \
  polymec_free(map); \
} \
\
static inline int map_name##_hash(map_name##_t* map, key_type key) \
{ \
  int h = map->hash(key); \
  h += ~(h << 9); \
  h ^= (((unsigned int) h) >> 14); \
  h += (h << 4); \
  h ^= (((unsigned int) h) >> 10); \
  return h; \
} \
\
static inline int map_name##_index(int bucket_count, int hash) \
{ \
  return hash & (bucket_count - 1); \
} \
\
static inline bool map_name##_keys_equal(map_name##_t* map, key_type key1, int hash1, key_type key2, int hash2) \
{ \
  return ((hash1 == hash2) && map->equals(key1, key2)); \
} \
\
static inline map_name##_entry_t* map_name##_get_entry(map_name##_t* map, key_type key) \
{ \
  int h = map_name##_hash(map, key); \
  int index = map_name##_index(map->bucket_count, h); \
  map_name##_entry_t* entry = map->buckets[index]; \
  while (entry != NULL) \
  { \
    if (map_name##_keys_equal(map, entry->key, entry->hash, key, h)) \
      return entry; \
    entry = entry->next; \
  } \
  return NULL; \
} \
\
static inline map_name##_value_t* map_name##_get(map_name##_t* map, key_type key) \
{ \
  map_name##_entry_t* entry = map_name##_get_entry(map, key); \
  if (entry != NULL) \
    return &(entry->value); \
  else \
    return NULL; \
} \
\
static inline bool map_name##_contains(map_name##_t* map, key_type key) \
{ \
  int h = map_name##_hash(map, key); \
  int index = map_name##_index(map->bucket_count, h); \
  map_name##_entry_t* entry = map->buckets[index]; \
  while (entry != NULL) \
  { \
    if (map_name##_keys_equal(map, entry->key, entry->hash, key, h)) \
      return true; \
    entry = entry->next; \
  } \
  return false; \
} \
\
static inline void map_name##_expand(map_name##_t* map) \
{ \
  if (map->size > (map->bucket_count * 3/4)) \
  { \
    int new_count = map->bucket_count * 2; \
    map_name##_entry_t** new_buckets = (map_name##_entry_t**)polymec_calloc(new_count, sizeof(map_name##_entry_t*)); \
    if (new_buckets == NULL) \
      return; \
    for (int i = 0; i < map->bucket_count; ++i) \
    { \
      map_name##_entry_t* entry = map->buckets[i]; \
      while (entry != NULL) \
      { \
        map_name##_entry_t* next = entry->next; \
        int index = map_name##_index(new_count, entry->hash); \
        entry->next = new_buckets[index]; \
        new_buckets[index] = entry; \
        entry = next; \
      } \
    } \
    polymec_free(map->buckets); \
    map->buckets = new_buckets; \
    map->bucket_count = new_count; \
  } \
} \
\
static inline void map_name##_insert_with_dtors(map_name##_t* map, key_type key, value_type value, map_name##_kv_dtor kv_dtor, map_name##_k_dtor k_dtor, map_name##_v_dtor v_dtor) \
{ \
  ASSERT((kv_dtor == NULL) || (k_dtor == NULL) || (v_dtor == NULL)); \
  int h = map_name##_hash(map, key); \
  int index = map_name##_index(map->bucket_count, h); \
  int depth = 0; \
  map_name##_entry_t** p = &(map->buckets[index]); \
  while (true) \
  { \
    map_name##_entry_t* current = *p; \
    if (current == NULL) \
    { \
      *p = (map_name##_entry_t*)polymec_malloc(sizeof(map_name##_entry_t)); \
      (*p)->key = key; \
      (*p)->hash = h; \
      (*p)->value = value; \
      (*p)->kv_dtor = kv_dtor; \
      (*p)->k_dtor = k_dtor; \
      (*p)->v_dtor = v_dtor; \
      (*p)->next = NULL; \
      map->size++; \
      map_name##_expand(map); \
      map->max_depth = (depth+1 > map->max_depth) ? depth+1 : map->max_depth; \
      return; \
    } \
    if (map_name##_keys_equal(map, current->key, current->hash, key, h)) \
    { \
      if (current->k_dtor != NULL) \
        current->k_dtor(current->key); \
      if (current->v_dtor != NULL) \
        current->v_dtor(current->value); \
      if (current->kv_dtor != NULL) \
        current->kv_dtor(current->key, current->value); \
      (*p)->key = key; \
      (*p)->kv_dtor = kv_dtor; \
      (*p)->k_dtor = k_dtor; \
      (*p)->v_dtor = v_dtor; \
      current->value = value; \
      return; \
    } \
    depth++; \
    p = &current->next; \
  } \
} \
\
static inline void map_name##_insert_with_kv_dtor(map_name##_t* map, key_type key, value_type value, map_name##_kv_dtor dtor) \
{ \
  map_name##_insert_with_dtors(map, key, value, dtor, NULL, NULL); \
} \
\
static inline void map_name##_insert_with_kv_dtors(map_name##_t* map, key_type key, value_type value, map_name##_k_dtor k_dtor, map_name##_v_dtor v_dtor) \
{ \
  map_name##_insert_with_dtors(map, key, value, NULL, k_dtor, v_dtor); \
} \
\
static inline void map_name##_insert_with_k_dtor(map_name##_t* map, key_type key, value_type value, map_name##_k_dtor dtor) \
{ \
  map_name##_insert_with_dtors(map, key, value, NULL, dtor, NULL); \
} \
\
static inline void map_name##_insert_with_v_dtor(map_name##_t* map, key_type key, value_type value, map_name##_v_dtor dtor) \
{ \
  map_name##_insert_with_dtors(map, key, value, NULL, NULL, dtor); \
} \
\
static inline void map_name##_insert(map_name##_t* map, key_type key, value_type value) \
{ \
  map_name##_insert_with_dtors(map, key, value, NULL, NULL, NULL); \
} \
\
static inline key_type map_name##_change_key(map_name##_t* map, key_type old_key, key_type new_key) \
{ \
  int h = map_name##_hash(map, old_key); \
  int index = map_name##_index(map->bucket_count, h); \
  map_name##_entry_t** p = &(map->buckets[index]); \
  map_name##_entry_t* current; \
  key_type key = old_key; \
  map_name##_kv_dtor kv_dtor = NULL; \
  map_name##_k_dtor k_dtor = NULL; \
  map_name##_v_dtor v_dtor = NULL; \
  while ((current = *p) != NULL) \
  { \
    if (map_name##_keys_equal(map, current->key, current->hash, old_key, h)) \
    { \
      *p = current->next; \
      key = current->key; \
      value_type value = current->value; \
      kv_dtor = current->kv_dtor; \
      k_dtor = current->k_dtor; \
      v_dtor = current->v_dtor; \
      polymec_free(current); \
      map->size--; \
      map_name##_insert_with_dtors(map, new_key, value, kv_dtor, k_dtor, v_dtor); \
      break; \
    }\
    else \
      p = &current->next; \
  } \
  return key; \
} \
\
static inline void map_name##_swap(map_name##_t* map, key_type key1, key_type key2) \
{ \
  map_name##_entry_t* e1 = map_name##_get_entry(map, key1); \
  map_name##_entry_t* e2 = map_name##_get_entry(map, key2); \
  ASSERT((e1 != NULL) && (e2 != NULL)); \
  ASSERT(e1->kv_dtor == e2->kv_dtor); \
  value_type temp_val = e1->value; \
  map_name##_v_dtor temp_v_dtor = e1->v_dtor; \
  e1->value = e2->value; e1->v_dtor = e2->v_dtor; \
  e2->value = temp_val; e2->v_dtor = temp_v_dtor; \
} \
\
static inline void map_name##_delete(map_name##_t* map, key_type key) \
{ \
  int h = map_name##_hash(map, key); \
  int index = map_name##_index(map->bucket_count, h); \
  map_name##_entry_t** p = &(map->buckets[index]); \
  map_name##_entry_t* current; \
  while ((current = *p) != NULL) \
  { \
    if (map_name##_keys_equal(map, current->key, current->hash, key, h)) \
    { \
      *p = current->next; \
      if (current->kv_dtor != NULL) \
        (*current->kv_dtor)(current->key, current->value); \
      else \
      { \
        if (current->k_dtor != NULL) \
          (*current->k_dtor)(current->key); \
        if (current->v_dtor != NULL) \
          (*current->v_dtor)(current->value); \
      } \
      polymec_free(current); \
      map->size--; \
    }\
    else \
      p = &current->next; \
  } \
} \
\
static inline bool map_name##_next(map_name##_t* map, int* pos, key_type* key, value_type* value) \
{ \
  int index = *pos / map->max_depth; \
  int depth = *pos % map->max_depth; \
  map_name##_entry_t* entry; \
  while ((index < map->bucket_count) && (map->buckets[index] == NULL)) index++; \
  if (index == map->bucket_count) \
    return false; \
  entry = map->buckets[index]; \
  for (int d = 0; d < depth; ++d) \
    entry = entry->next; \
  *key = entry->key; \
  *value = entry->value; \
  if (entry->next != NULL) \
    (*pos)++; \
  else \
    *pos = (index+1) * map->max_depth; \
  return true; \
} \
\
static inline map_name##_t* map_name##_clone(map_name##_t* map, \
                                             map_name##_key_t (*key_clone_func)(map_name##_key_t), \
                                             map_name##_value_t (*val_clone_func)(map_name##_value_t), \
                                             map_name##_k_dtor k_dtor, \
                                             map_name##_v_dtor v_dtor) \
{ \
  map_name##_t* clone = map_name##_new_with_capacity(map->bucket_count); \
  int pos = 0; \
  map_name##_key_t key; \
  map_name##_value_t value; \
  while (map_name##_next(map, &pos, &key, &value)) \
  { \
    map_name##_key_t key_clone; \
    map_name##_value_t val_clone; \
    if (key_clone_func != NULL) \
      key_clone = key_clone_func(key); \
    else \
      key_clone = key; \
    if (val_clone_func != NULL) \
      val_clone = val_clone_func(value); \
    else \
      val_clone = value; \
    map_name##_insert_with_kv_dtors(clone, key_clone, val_clone, k_dtor, v_dtor); \
  } \
  return clone; \
} \
\
static inline bool map_name##_empty(map_name##_t* map) \
{ \
  return (map->size == 0); \
} \
\

///@}

// Define some unordered maps.

/// \class int_int_unordered_map
/// A mapping of integers to integers.
DEFINE_UNORDERED_MAP(int_int_unordered_map, int, int, int_hash, int_equals)

/// \class int_index_unordered_map
/// A mapping of integers to 64-bit indices.
DEFINE_UNORDERED_MAP(int_index_unordered_map, int, index_t, int_hash, int_equals)

/// \class int_real_unordered_map
/// A mapping of integers to real numbers.
DEFINE_UNORDERED_MAP(int_real_unordered_map, int, real_t, int_hash, int_equals)

/// \class int_ptr_unordered_map
/// A mapping of integers to pointers.
DEFINE_UNORDERED_MAP(int_ptr_unordered_map, int, void*, int_hash, int_equals)

/// \class int_pair_int_unordered_map
/// A mapping of integer pairs (int[2]) to integers.
DEFINE_UNORDERED_MAP(int_pair_int_unordered_map, int*, int, int_pair_hash, int_pair_equals)

/// \class int_tuple_int_unordered_map
/// A mapping of integer tuples (int*) to integers.
DEFINE_UNORDERED_MAP(int_tuple_int_unordered_map, int*, int, int_tuple_hash, int_tuple_equals)

/// \class index_index_unordered_map
/// A mapping of 64-bit indices to 64-bit indices.
DEFINE_UNORDERED_MAP(index_index_unordered_map, index_t, index_t, index_hash, index_equals)

/// \class index_int_unordered_map
/// A mapping of 64-bit indices to integers.
DEFINE_UNORDERED_MAP(index_int_unordered_map, index_t, int, index_hash, index_equals)

/// \class index_real_unordered_map
/// A mapping of 64-bit indices to real numbers.
DEFINE_UNORDERED_MAP(index_real_unordered_map, index_t, real_t, index_hash, index_equals)

/// \class index_pair_int_unordered_map
/// A mapping of 64-bit index pairs (index_t[2]) to integers.
DEFINE_UNORDERED_MAP(index_pair_int_unordered_map, index_t*, int, index_pair_hash, index_pair_equals)

/// \class index_pair_real_unordered_map
/// A mapping of 64-bit index pairs (index_t[2]) to real numbers.
DEFINE_UNORDERED_MAP(index_pair_real_unordered_map, index_t*, real_t, index_pair_hash, index_pair_equals)

/// \class index_ptr_unordered_map
/// A mapping of 64-bit indices to real numbers.
DEFINE_UNORDERED_MAP(index_ptr_unordered_map, index_t, void*, index_hash, index_equals)

/// \class string_int_unordered_map
/// A mapping of strings (char*) to integers.
DEFINE_UNORDERED_MAP(string_int_unordered_map, char*, int, string_hash, string_equals)

/// \class string_string_unordered_map
/// A mapping of strings (char*) to strings.
DEFINE_UNORDERED_MAP(string_string_unordered_map, char*, char*, string_hash, string_equals)

/// \class string_ptr_unordered_map
/// A mapping of strings (char*) to pointers.
DEFINE_UNORDERED_MAP(string_ptr_unordered_map, char*, void*, string_hash, string_equals)

/// \class string_ptr_unordered_map
/// A mapping of pointers to pointers.
DEFINE_UNORDERED_MAP(ptr_ptr_unordered_map, void*, void*, ptr_hash, ptr_equals)

#endif
