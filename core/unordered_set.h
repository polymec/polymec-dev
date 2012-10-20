#ifndef ARBI_UNORDERED_SET_H
#define ARBI_UNORDERED_SET_H

#include "core/hash_map.h"
#include "core/hash_functions.h"
#include "core/comparators.h"

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
// void x_set_foreach(x_set_t* set, x_set_visitor visit, void*) - Executes visit on each set element.
// void x_set_copy_out(x_set_t* set, x* array) - Copies all elements from the set to the given array.

#define DEFINE_UNORDERED_SET(set_name, element, hash_func, equals_func) \
DEFINE_HASH_MAP(set_name##_hash_map, element, bool, hash_func, equals_func) \
typedef element set_name##_element_t; \
typedef int (*set_name##_hash_func)(element); \
typedef struct \
{ \
  set_name##_hash_map_t* map; \
  int size; \
} set_name##_t; \
\
static inline set_name##_t* set_name##_new() \
{ \
  set_name##_t* set = malloc(sizeof(set_name##_t)); \
  set->map = set_name##_hash_map_new(); \
  set->size = 0; \
  return set; \
} \
\
static inline void set_name##_clear(set_name##_t* set) \
{ \
  set_name##_hash_map_clear(set->map); \
  set->size = 0; \
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
  return set_name##_hash_map_contains(set->map, datum); \
} \
\
static inline void set_name##_insert(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_hash_map_set(set->map, datum, true); \
  set->size = set->map->size; \
} \
\
static inline void set_name##_delete(set_name##_t* set, set_name##_element_t datum) \
{ \
  set_name##_hash_map_delete(set->map, datum); \
  set->size = set->map->size; \
} \
\
typedef struct \
{ \
  int index; \
  element* array; \
} set_name##_hash_copy_out_t; \
static inline void set_name##_hash_copy_out(element key, bool value, void* arg) \
{ \
  set_name##_hash_copy_out_t* a = (set_name##_hash_copy_out_t*)arg; \
  a->array[a->index] = key; \
  a->index++; \
} \
\
static inline void set_name##_copy_out(set_name##_t* set, element* array) \
{ \
  set_name##_hash_copy_out_t arg = {.index = 0, .array = array}; \
  set_name##_hash_map_foreach(set->map, set_name##_hash_copy_out, (void*)&arg); \
} \
\

// Define some ordered_sets.
DEFINE_UNORDERED_SET(int_unordered_set, int, int_hash, int_equals)
DEFINE_UNORDERED_SET(string_unordered_set, char*, string_hash, string_equals)

#endif
