// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ORDERED_SET_H
#define POLYMEC_ORDERED_SET_H

#include "core/avl_tree.h"

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
  set_name##_t* set = (set_name##_t*)polymec_malloc(sizeof(set_name##_t)); \
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
  polymec_free(set); \
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

/// \addtogroup core core
///@{

/// \def DEFINE_ORDERED_SET(set_name, element, comparator)
/// Defines an ordered set for the given element type. The following
/// interface is defined for a set with map_name `x_set`.
/// * `x_set_t* x_set_new()` - Creates a new empty ordered set.
/// * `void x_set_free(x_set_t* set)` - Destroys the set.
/// * `void x_set_clear(x_set_t* set)` - Empties the set.
/// * `x_set_node_t* x_set_find(x_set_t* set, x datum)` - Finds the node containing the datum.
/// * `bool x_set_contains(x_set_t* set, x datum)` - Returns true if the set contains the datum, false otherwise.
/// * `void x_set_insert(x_set_t* set, x datum)` - Inserts a datum into the set.
/// * `void x_set_delete(x_set_t* set, x datum)` - Deletes the datum from the set.
///
/// Member data for a set `set`:
/// * set->size - The number of element in the set
///
/// \param set_name The name of the ordered set.
/// \param element The data type stored by the set.
/// \param comparator A comparator function that takes elements a, b, and
///                   returns -1 if a < b, 0 if a == b, and 1 if a > b.
#define DEFINE_ORDERED_SET(set_name, element, comparator) \
DEFINE_AVL_TREE(set_name##_avl_tree, element, comparator) \
DEFINE_ORDERED_SET_USING_TREE(set_name, set_name##_avl_tree)

///@}

// Define some ordered_sets.
DEFINE_ORDERED_SET_USING_AVL_TREE(int_ordered_set, int_avl_tree)
DEFINE_ORDERED_SET_USING_AVL_TREE(real_ordered_set, real_avl_tree)
DEFINE_ORDERED_SET_USING_AVL_TREE(string_ordered_set, string_avl_tree)

#endif
