#ifndef ARBI_ORDERED_MAP_H
#define ARBI_ORDERED_MAP_H

#include "core/avl_tree.h"
#include "core/key_value.h"

// An ordered map is a container that associates keys with values.
// One defines an ordered map and node types using
// DEFINE_ORDERED_MAP(map_name, key_type, value_type, comparator, destructor)
// or 
// DEFINE_ORDERED_MAP_USING_AVL_TREE(map_name, tree_type),
// which uses tree_type for the underlying implementation.

// Interface for a type x_map_t (with node type x_map_node_t 
// and datum x) defined with 
// DEFINE_ORDERED_MAP(x_map, x, x_comparator, x_destructor):
// 
// x_map_t* x_map_new() - Creates a new empty ordered set.
// void x_map_free(x_map_t* set) - Destroys the set.
// void x_map_clear(x_map_t* set) - Empties the set.
// x_map_node_t* x_map_find(x_map_t* set, x datum) - Finds the node containing the datum.
// bool x_map_contains(x_map_t* set, x datum) - Returns true if the set contains the datum, false otherwise.
// void x_map_insert(x_map_t* set, x datum) - Inserts a datum into the set.
// void x_map_delete(x_map_t* set, x datum) - Deletes the datum from the set.
// void x_map_foreach(x_map_t* node, x_map_visitor visit, void*) - Executes visit on each set element.

#define DEFINE_ORDERED_MAP_USING_TREE(map_name, tree_name) \
typedef struct map_name##_node_t map_name##_node_t; \
typedef void (*map_name##_node_visitor)(map_name##_node_t*, void*); \
typedef tree_name##_node_t map_name##_node_t; \
\
typedef struct map_name##_t map_name##_t; \
typedef tree_name##_destructor map_name##_destructor; \
struct map_name##_t \
{ \
  tree_name##_t tree; \
  int size; \
}; \
\
static inline map_name##_t* map_name##_new() \
{ \
  map_name##_t* set = malloc(sizeof(map_name##_t)); \
  set->tree = tree_name##_new(); \
  set->size = 0; \
  return set; \
} \
\
static inline void map_name##_clear(map_name##_t* set) \
{ \
  tree_name##_clear(set->tree); \
  set->size = 0; \
} \
\
static inline void map_name##_free(map_name##_t* set) \
{ \
  tree_name##_free(set->tree); \
  free(set); \
} \
\
static inline map_name##_node_t* map_name##_find(map_name##_t* set, element datum) \
{ \
  return tree_name##_find(set->tree, datum); \
} \
\
static inline bool map_name##_contains(map_name##_t* set, element datum) \
{ \
  return (map_name##_find(set, datum) != NULL); \
} \
\
static inline void map_name##_insert(map_name##_t* set, element datum) \
{ \
  tree_name##_insert(set->tree, datum); \
  set->size = tree_name##_size(set->tree); \
} \
\
static inline void map_name##_delete(map_name##_t* set, element datum) \
{ \
  map_name##_node_t* node = map_name##_find(set, datum); \
  if (node != NULL) \
  { \
    tree_name##_delete(set->tree, node); \
    set->size = tree_name##_size(set->tree); \
  } \
} \
\
static inline void map_name##_foreach(map_name##_t* set, map_name##_visitor visit, void* arg) \
{ \
  tree_name##_node_visit(set->tree->root, visit, arg); \
} \
\

#define DEFINE_ORDERED_MAP(map_name, key_type, value_type, comparator, destructor) \
DEFINE_KEY_VALUE(map_name##_key_value, key_type, value_type) \
DEFINE_AVL_TREE(map_name##_avl_tree, map_name##_key_value, comparator, destructor) \
DEFINE_ORDERED_SET_USING_TREE(map_name, map_name##_avl_tree)

#endif
