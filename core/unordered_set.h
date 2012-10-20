#ifndef ARBI_UNORDERED_SET_H
#define ARBI_UNORDERED_SET_H

#include "core/hash_map.h"

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
// void x_set_foreach(x_set_t* node, x_set_visitor visit, void*) - Executes visit on each set element.

#define DEFINE_UNORDERED_SET(set_name, element, hash_func, equals_func) \
DEFINE_HASH_MAP(set_name##_hash_map, element, bool, hash_func, equals_func) \
typedef element set_name##_element_t; \
typedef int (*set_name##_hash_func)(element); \
typedef struct \
{ \
  set_name##_hash_map map; \
} set_name##_t; \
\
static inline set_name##_t* set_name##_new() \
{ \
  set_name##_t* set = malloc(sizeof(set_name##_t)); \
  set->map = set_name##_hash_map_new(); \
  return set; \
} \
\
static inline void set_name##_clear(set_name##_t* set) \
{ \
  set_name##_hash_map_clear(set->map); \
} \
\
static inline void set_name##_free(set_name##_t* set) \
{ \
  set_name##_hash_map_free(set->map); \
  free(set); \
} \
\
static inline bool set_name##_contains(set_name##_t* set, set_name##_element_t datum) \
{ \
  return set_name##_hash_map_contains(set, datum); \
} \
\
static inline void set_name##_insert(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_hash_map_set(set->map, datum, true); \
} \
\
static inline void set_name##_delete(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_hash_map_delete(set->map, datum); \
} \
\

// Define some ordered_sets.
DEFINE_UNORDERED_SET(int_unordered_set, int, int_hash, int_equals)
DEFINE_UNORDERED_SET(double_unordered_set, double, double_hash, double_equals)
DEFINE_UNORDERED_SET(string_unordered_set, char*, string_hash, double_hash)

#endif
