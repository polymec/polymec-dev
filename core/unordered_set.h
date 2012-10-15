#ifndef ARBI_UNORDERED_SET_H
#define ARBI_UNORDERED_SET_H

// An unordered set is a container that stores unique values using a hash table.
// One defines an unordered set and node types using
// DEFINE_UNORDERED_SET(set_name, element, hash_function, destructor)
//
// Interface for a type x_set_t (with node type x_set_node_t 
// and datum x) defined with 
// DEFINE_UNORDERED_SET(x_set, x, x_comparator, x_destructor):
// 
// x_set_t* x_set_new() - Creates a new empty ordered set.
// void x_set_free(x_set_t* set) - Destroys the set.
// void x_set_clear(x_set_t* set) - Empties the set.
// x_set_node_t* x_set_find(x_set_t* set, x datum) - Finds the node containing the datum.
// bool x_set_contains(x_set_t* set, x datum) - Returns true if the set contains the datum, false otherwise.
// void x_set_insert(x_set_t* set, x datum) - Inserts a datum into the set.
// void x_set_delete(x_set_t* set, x datum) - Deletes the datum from the set.
// void x_set_foreach(x_set_t* node, x_set_visitor visit, void*) - Executes visit on each set element.

#define DEFINE_UNORDERED_SET(set_name, element, hash_function, destructor) \
typedef int (*set_name##_hash_function)(element); \
typedef struct set_name##_t set_name##_t; \
typedef void (*set_name##_destructor)(element); \
typedef struct \
{ \
  int hash; \
  element value; \
} set_name##_entry_t; \
struct set_name##_t \
{ \
  set_name##_entry_t* entries; \
  set_name##_hash_function hash; \
  int size; \
}; \
\
static inline set_name##_t* set_name##_new() \
{ \
  set_name##_t* set = malloc(sizeof(set_name##_t)); \
  set->entries = NULL; \
  set->hash = hash_function; \
  set->size = 0; \
  return set; \
} \
\
static inline void set_name##_clear(set_name##_t* set) \
{ \
  if (set->entries != NULL) \
  { \
    free(set->entries); \
    set->entries = NULL; \
  } \
  set->size = 0; \
} \
\
static inline void set_name##_free(set_name##_t* set) \
{ \
  set_name##_clear(set); \
  free(set); \
} \
\
static inline set_name##_entry_t* set_name##_find(set_name##_t* set, set_name##_element_t datum) \
{ \
  int key_hash = set->hash(datum); \
  for (int i = 0; i < set->size; ++i) \
  { \
  } \
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
static inline void set_name##_foreach(set_name##_t* set, set_name##_visitor visit, void* arg) \
{ \
  tree_name##_node_visit(set->tree->root, visit, arg); \
} \
\

#define DEFINE_ORDERED_SET(set_name, element, comparator, destructor) \
DEFINE_AVL_TREE(set_name##_avl_tree, element, comparator, destructor) \
DEFINE_ORDERED_SET_USING_TREE(set_name, set_name##_avl_tree)

// Define some ordered_sets.
DEFINE_ORDERED_SET_USING_AVL_TREE(int_ordered_set, int_avl_tree)
DEFINE_ORDERED_SET_USING_AVL_TREE(double_ordered_set, double_avl_tree)
DEFINE_ORDERED_SET_USING_AVL_TREE(string_ordered_set, string_avl_tree)

#endif
