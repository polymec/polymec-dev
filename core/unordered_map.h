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

#ifndef POLYMEC_UNORDERED_MAP_H
#define POLYMEC_UNORDERED_MAP_H

#include "core/polymec.h"
#include "core/hash_functions.h"
#include "core/comparators.h"
#include "core/tuple.h"

// An unordered map is a map that associates a key with a single value, 
// using a hash to store them. One defines an unordered map using
// DEFINE_UNORDERED_MAP(map_name, key_type, value_type, hash_func, key_equals).

// Interface for a type x_map_t (with key type K and 
// value type V) defined with 
// DEFINE_UNORDERED_MAP(x_map, K, V, x_hash, key_equals):
// 
// x_map_t* x_map_new() - Creates a new empty unordered map.
// x_map_t* x_map_new_with_capacity(int N) - Creates a new empty hash map with 
//                                           initial capacity N.
// x_map_t* x_map_copy(x_map_t* m) - Creates a (shallow) copy of m.
// void x_map_free(x_map_t* map) - Destroys the map.
// void x_map_clear(x_map_t* map) - Empties the map.
// x_map_value_t* x_map_get(x_map_t* map, x_map_key_t key) - Returns the value for the key, or NULL.
// bool x_map_contains(x_map_t* map, x_map_key_t key) - Returns true if the map contains the key, false if not.
// void x_map_insert(x_map_t* map, x_map_key_t key, x_map_value_t value) - Sets the value for the given key.
// void x_map_insert_with_kv_dtor(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_kv_dtor) - Sets the value for the given key and associates a destructor for the key/value pair.
// void x_map_insert_with_k_dtor(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_k_dtor) - Sets the value for the given key and associates a destructor for the key.
// void x_map_insert_with_v_dtor(x_map_t* map, x_map_key_t key, x_map_value_t value, x_map_v_dtor) - Sets the value for the given key and associates a destructor for the value.
// key_type x_map_change_key(x_map_t* map, x_map_key_t old_key, x_map_key_t new_key) - Renames old_key to new_key, overwriting new_key if it exists. Returns old key.
// void x_map_delete(x_map_t* map, x_map_key_t key) - Deletes the value for the given key.
// bool x_map_next(x_map_t* map, int* pos, x_map_key_t* key, x_map_value_t* value) - Allows the traversal of the maps keys and values.
// bool x_map_empty(x_map_t* map) - Returns true if empty, false otherwise.

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
  map->buckets = (map_name##_entry_t**)polymec_malloc(map->bucket_count * sizeof(map_name##_entry_t*)); \
  memset(map->buckets, 0, map->bucket_count * sizeof(map_name##_entry_t*)); \
  ASSERT(map->buckets != NULL); \
  map->size = 0; \
  map->hash = hash_func; \
  map->equals = equals_func; \
  map->max_depth = 1; \
  return map; \
} \
\
static inline map_name##_t* map_name##_new() \
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
      else if (entry->k_dtor != NULL) \
        (*entry->k_dtor)(entry->key); \
      else if (entry->v_dtor != NULL) \
        (*entry->v_dtor)(entry->value); \
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
static inline map_name##_value_t* map_name##_get(map_name##_t* map, key_type key) \
{ \
  int h = map_name##_hash(map, key); \
  int index = map_name##_index(map->bucket_count, h); \
  map_name##_entry_t* entry = map->buckets[index]; \
  while (entry != NULL) \
  { \
    if (map_name##_keys_equal(map, entry->key, entry->hash, key, h)) \
      return &(entry->value); \
    entry = entry->next; \
  } \
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
    map_name##_entry_t** new_buckets = (map_name##_entry_t**)polymec_malloc(new_count * sizeof(map_name##_entry_t*)); \
    memset(new_buckets, 0, new_count * sizeof(map_name##_entry_t*)); \
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
  value_type value; \
  map_name##_kv_dtor kv_dtor = NULL; \
  map_name##_k_dtor k_dtor = NULL; \
  map_name##_v_dtor v_dtor = NULL; \
  bool found = false; \
  while ((current = *p) != NULL) \
  { \
    if (map_name##_keys_equal(map, current->key, current->hash, old_key, h)) \
    { \
      *p = current->next; \
      key = current->key; \
      value = current->value; \
      kv_dtor = current->kv_dtor; \
      k_dtor = current->k_dtor; \
      v_dtor = current->v_dtor; \
      polymec_free(current); \
      map->size--; \
      found = true; \
    }\
    else \
      p = &current->next; \
  } \
  if (found) \
    map_name##_insert_with_dtors(map, new_key, value, kv_dtor, k_dtor, v_dtor); \
  return key; \
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
      else if (current->k_dtor != NULL) \
        (*current->k_dtor)(current->key); \
      else if (current->v_dtor != NULL) \
        (*current->v_dtor)(current->value); \
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
static inline map_name##_t* map_name##_copy(map_name##_t* map) \
{ \
  map_name##_t* copy = map_name##_new_with_capacity(map->bucket_count); \
  int pos = 0; \
  map_name##_key_t key; \
  map_name##_value_t value; \
  while (map_name##_next(map, &pos, &key, &value)) \
    map_name##_insert(copy, key, value); \
  return copy; \
} \
\
static inline bool map_name##_empty(map_name##_t* map) \
{ \
  return (map->size == 0); \
} \
\

// Define some unordered maps.
DEFINE_UNORDERED_MAP(int_int_unordered_map, int, int, int_hash, int_equals)
DEFINE_UNORDERED_MAP(int_index_unordered_map, int, index_t, int_hash, int_equals)
DEFINE_UNORDERED_MAP(int_ptr_unordered_map, int, void*, int_hash, int_equals)
DEFINE_UNORDERED_MAP(int_pair_int_unordered_map, int*, int, int_pair_hash, int_pair_equals)
DEFINE_UNORDERED_MAP(int_tuple_int_unordered_map, int*, int, int_tuple_hash, int_tuple_equals)
DEFINE_UNORDERED_MAP(index_index_unordered_map, index_t, index_t, index_hash, index_equals)
DEFINE_UNORDERED_MAP(index_int_unordered_map, index_t, int, index_hash, index_equals)
DEFINE_UNORDERED_MAP(index_ptr_unordered_map, index_t, void*, index_hash, index_equals)
DEFINE_UNORDERED_MAP(string_int_unordered_map, char*, int, string_hash, string_equals)
DEFINE_UNORDERED_MAP(string_string_unordered_map, char*, char*, string_hash, string_equals)
DEFINE_UNORDERED_MAP(string_ptr_unordered_map, char*, void*, string_hash, string_equals)
DEFINE_UNORDERED_MAP(ptr_ptr_unordered_map, void*, void*, pointer_hash, pointer_equals)

#endif
