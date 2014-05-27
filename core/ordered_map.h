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

#ifndef POLYMEC_ORDERED_MAP_H
#define POLYMEC_ORDERED_MAP_H

#include "core/avl_tree.h"
#include "core/key_value.h"

// An ordered map is a container that associates keys with values.
// One defines an ordered map and node types using
// DEFINE_ORDERED_MAP(map_name, key_type, value_type, comparator).

// Interface for a type x_map_t (with node type x_map_node_t 
// and datum x) defined with 
// DEFINE_ORDERED_MAP(x_map, x, x_comparator):
// 
// x_map_t* x_map_new() - Creates a new empty ordered map.
// void x_map_free(x_map_t* map) - Destroys the map.
// void x_map_clear(x_map_t* map) - Empties the map.
// x_map_node_t* x_map_find(x_map_t* map, x_map_key_t key) - Finds the node whose key matches the given one.
// x_map_value_t x_map_value(x_map_t* map, x_map_key_t key) - Finds the value for the given key.
// void x_map_insert(x_map_t* map, x datum) - Inserts a datum into the map.
// void x_map_delete(x_map_t* map, x datum) - Deletes the datum from the map.

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
  map_name##_t* map = polymec_malloc(sizeof(map_name##_t)); \
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

#endif
