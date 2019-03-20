// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ORDERED_MAP_H
#define POLYMEC_ORDERED_MAP_H

#include "core/avl_tree.h"
#include "core/key_value.h"

/// \addtogroup core core
///@{

/// \def DEFINE_ORDERED_MAP(map_name, key_type, value_type, key_comparator)
/// Defines an ordered map for the given key and value types. The following
/// interface is defined for a map with map_name `x_map`.
/// * `x_map_t* x_map_new()` - Creates a new empty ordered map.
/// * `void x_map_free(x_map_t* map)` - Destroys the map.
/// * `void x_map_clear(x_map_t* map)` - Empties the map.
/// * `x_map_node_t* x_map_find(x_map_t* map, x_map_key_t key)` - Finds the node whose key matches the given one.
/// * `x_map_value_t x_map_value(x_map_t* map, x_map_key_t key)` - Finds the value for the given key.
/// * `void x_map_insert(x_map_t* map, x datum)` - Inserts a datum into the map.
/// * `void x_map_delete(x_map_t* map, x datum)` - Deletes the datum from the map.
///
/// Member data for an ordered map `map`:
/// * `map->size` - The number of key-value pairs in the map.
///
/// \param map_name The name of the unordered map.
/// \param key_type The data type used as a key in the map.
/// \param value_type The data type used as a value in the map.
/// \param key_comparator A comparator function that accepts two key arguments a, b, and
///                       returns -1 if a < b, 0 if a == b, and 1 if a > b.
#define DEFINE_ORDERED_MAP(map_name, key_type, value_type, key_comparator) \
DEFINE_KEY_VALUE(map_name##_key_value, key_type, value_type) \
static inline int map_name##_key_value_cmp(map_name##_key_value_t x, map_name##_key_value_t y) \
{ \
  return key_comparator(x.key, y.key); \
} \
DEFINE_AVL_TREE(map_name##_avl_tree, map_name##_key_value_t, map_name##_key_value_cmp) \
typedef map_name##_avl_tree##_node_t map_name##_node_t; \
typedef key_type map_name##_key_t; \
typedef value_type map_name##_value_t; \
typedef void (*map_name##_visitor)(map_name##_node_t*, void*); \
typedef struct map_name##_t map_name##_t; \
struct map_name##_t \
{ \
  map_name##_avl_tree_t* tree; \
  int size; \
}; \
\
static inline map_name##_t* map_name##_new() \
{ \
  map_name##_t* map = (map_name##_t*)polymec_malloc(sizeof(map_name##_t)); \
  map->tree = map_name##_avl_tree_new(); \
  map->size = 0; \
  return map; \
} \
\
static inline void map_name##_clear(map_name##_t* map) \
{ \
  map_name##_avl_tree_clear(map->tree); \
  map->size = 0; \
} \
\
static inline void map_name##_free(map_name##_t* map) \
{ \
  map_name##_avl_tree_free(map->tree); \
  polymec_free(map); \
} \
\
static inline map_name##_node_t* map_name##_find(map_name##_t* map, key_type key) \
{ \
  map_name##_key_value_t kv = {.key = key}; \
  return map_name##_avl_tree_find(map->tree, kv); \
} \
\
static inline map_name##_value_t map_name##_value(map_name##_t* map, key_type key) \
{ \
  map_name##_node_t* node = map_name##_find(map, key); \
  map_name##_key_value_t kv; \
  if (node != NULL) \
  { \
    kv.key = key; \
    kv.value = node->value.value; \
  } \
  return kv.value; \
} \
\
static inline void map_name##_insert(map_name##_t* map, key_type key, value_type value) \
{ \
  map_name##_key_value_t kv = {.key = key, .value = value}; \
  map_name##_avl_tree_insert(map->tree, kv); \
  map->size = map_name##_avl_tree_size(map->tree); \
} \
\
static inline void map_name##_delete(map_name##_t* map, key_type key) \
{ \
  map_name##_node_t* node = map_name##_find(map, key); \
  if (node != NULL) \
  { \
    map_name##_avl_tree_delete(map->tree, node); \
    map->size = map_name##_avl_tree_size(map->tree); \
  } \
} \
\

///@}

#endif
