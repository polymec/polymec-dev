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

#ifndef POLYMEC_ORDERED_SET_H
#define POLYMEC_ORDERED_SET_H

#include "core/avl_tree.h"

// An ordered set is a container that contains unique, ordered values.
// One defines an ordered set and node types using
// DEFINE_ORDERED_SET(set_name, element, comparator)
// or 
// DEFINE_ORDERED_SET_USING_AVL_TREE(set_name, tree_type),
// which uses tree_type for the underlying implementation.

// Interface for a type x_set_t (with node type x_set_node_t 
// and datum x) defined with 
// DEFINE_ORDERED_SET(x_set, x, x_comparator):
// 
// x_set_t* x_set_new() - Creates a new empty ordered set.
// void x_set_free(x_set_t* set) - Destroys the set.
// void x_set_clear(x_set_t* set) - Empties the set.
// x_set_node_t* x_set_find(x_set_t* set, x datum) - Finds the node containing the datum.
// bool x_set_contains(x_set_t* set, x datum) - Returns true if the set contains the datum, false otherwise.
// void x_set_insert(x_set_t* set, x datum) - Inserts a datum into the set.
// void x_set_delete(x_set_t* set, x datum) - Deletes the datum from the set.

#define DEFINE_ORDERED_SET_USING_AVL_TREE(set_name, tree_name) \
typedef tree_name##_node_t set_name##_node_t; \
typedef tree_name##_element_t set_name##_element_t; \
typedef void (*set_name##_visitor)(set_name##_node_t*, void*); \
\
typedef struct set_name##_t set_name##_t; \
struct set_name##_t \
{ \
  tree_name##_t* tree; \
  int size; \
}; \
\
static inline set_name##_t* set_name##_new() \
{ \
  set_name##_t* set = malloc(sizeof(set_name##_t)); \
  set->tree = tree_name##_new(); \
  set->size = 0; \
  return set; \
} \
\
static inline void set_name##_clear(set_name##_t* set) \
{ \
  tree_name##_clear(set->tree); \
  set->size = 0; \
} \
\
static inline void set_name##_free(set_name##_t* set) \
{ \
  tree_name##_free(set->tree); \
  free(set); \
} \
\
static inline set_name##_node_t* set_name##_find(set_name##_t* set, set_name##_element_t datum) \
{ \
  return tree_name##_find(set->tree, datum); \
} \
\
static inline bool set_name##_contains(set_name##_t* set, set_name##_element_t datum) \
{ \
  return (set_name##_find(set, datum) != NULL); \
} \
\
static inline void set_name##_insert(set_name##_t* set, set_name##_element_t datum) \
{ \
  tree_name##_insert(set->tree, datum); \
  set->size = tree_name##_size(set->tree); \
} \
\
static inline void set_name##_delete(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_node_t* node = set_name##_find(set, datum); \
  if (node != NULL) \
  { \
    tree_name##_delete(set->tree, node); \
    set->size = tree_name##_size(set->tree); \
  } \
} \
\

#define DEFINE_ORDERED_SET(set_name, element, comparator) \
DEFINE_AVL_TREE(set_name##_avl_tree, element, comparator) \
DEFINE_ORDERED_SET_USING_TREE(set_name, set_name##_avl_tree)

// Define some ordered_sets.
DEFINE_ORDERED_SET_USING_AVL_TREE(int_ordered_set, int_avl_tree)
DEFINE_ORDERED_SET_USING_AVL_TREE(real_ordered_set, real_avl_tree)
DEFINE_ORDERED_SET_USING_AVL_TREE(string_ordered_set, string_avl_tree)

#endif
