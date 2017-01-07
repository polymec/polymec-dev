// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ARRAY_H
#define POLYMEC_ARRAY_H

#include "core/polymec.h"

// An array is a dynamically-resizable chunk of contiguous memory that 
// stores a particular data type. It serves the same purpose as a C++ vector.
// One defines an array using
// DEFINE_ARRAY(array_name, element)
//
// Interface for a type x_array_t (with datum x) defined with 
// DEFINE_ARRAY(x_array, x):
//
// x_array_t* x_array_new() - Creates a new, empty array on the heap.
// x_array_t* x_array_new_with_capacity(int capacity) - Creates a new, empty array on the heap with the given capacity.
// x_array_t empty_x_array() - Creates a new, empty array on the stack.
// void x_array_free(x_array_t* array) - Destroys the (heap-allocated) array.
// x* x_array_find(x_array_t* array, x value, cmp_func comparator) - Performs a linear search within the array, returning the pointer to the found item or NULL if not found.
// void x_array_append(x_array_t* array, x value) - Appends an x to the end of the array.
// void x_array_append_with_dtor(x_array_t* array, x value, destructor dtor) - Appends an x to the end of the array, using dtor to destroy when finished.
// void x_array_insert(x_array_t* array, size_t i, x value) - Inserts an x at position i within the array, resizing as necessary.
// void x_array_insert_with_dtor(x_array_t* array, size_t i, x value, destructor dtor) - Inserts an x at position i within the array, using dtor to destroy when finished.
// void x_array_assign(x_array_t* array, size_t i, x value) - Assigns an x to position i within the array.
// void x_array_assign_with_dtor(x_array_t* array, size_t i, x value, destructor dtor) - Assigns an x to position i within the array, using dtor to destroy when finished.
// bool x_array_empty(x_array_t* array) - Returns true if empty, false otherwise.
// void x_array_clear(x_array_t* array) - Clears the given array, making it empty.
// void x_array_resize(x_array_t* array, size_t new_size) - Resizes the array, keeping data intact if possible.
// void x_array_reserve(x_array_t* array, size_t new_capacity) - Reserves storage for the given capacity within the array. No effect if the array already has sufficient storage.
// bool x_array_next(x_array_t* array, int* pos, x* element) - Allows traversal over the items in the array.
// void x_array_release_data_and_free(x_array_t* array) - Releases control of the array data, allowing another entity to assume responsibility. Also destroys this object. USE CAREFULLY!!!
// 
// Member data for an array a:
// 
// a.data - A pointer to an array containing the actual data.
// a.dtors - A pointer to an array of destructors for data elements.
// a.size - Number of elements in the array.
// a.capacity - Number of elements in the array.

#define DEFINE_ARRAY(array_name, element) \
typedef void (*array_name##_dtor)(element); \
typedef struct \
{ \
  element* data; \
  array_name##_dtor* dtors; \
  size_t size; \
  size_t capacity; \
} array_name##_t; \
\
typedef int (*array_name##_comparator)(element, element); \
\
static inline void array_name##_reserve(array_name##_t* array, size_t new_capacity) \
{ \
  if (new_capacity > array->capacity) \
  { \
    array->data = (element*)polymec_realloc(array->data, sizeof(element) * new_capacity); \
    array->capacity = new_capacity; \
  } \
} \
\
static inline array_name##_t* array_name##_new_with_capacity(size_t capacity) \
{ \
  array_name##_t* array = (array_name##_t*)polymec_malloc(sizeof(array_name##_t)); \
  array->data = NULL; \
  array->dtors = NULL; \
  array->size = 0; \
  array->capacity = 0; \
  array_name##_reserve(array, capacity); \
  return array; \
} \
\
static inline array_name##_t* array_name##_new() \
{ \
  return array_name##_new_with_capacity(16); \
} \
\
static inline array_name##_t empty_##array_name() \
{ \
  static array_name##_t empty = {NULL, NULL, 0, 0}; \
  return empty; \
} \
\
static inline void array_name##_resize(array_name##_t* array, size_t new_size) \
{ \
  if (new_size > array->capacity) \
  { \
    while (array->capacity < new_size) \
      array_name##_reserve(array, MAX(16, 2 * array->capacity)); \
  } \
  else if (new_size < array->size) \
  { \
    if (array->dtors != NULL) \
    { \
      for (size_t i = new_size; i < array->size; ++i) \
      { \
        if (array->dtors[i] != NULL) \
        { \
          array->dtors[i](array->data[i]); \
          array->dtors[i] = NULL; \
        } \
      } \
    } \
  } \
  array->size = new_size; \
} \
\
static inline void array_name##_clear(array_name##_t* array) \
{ \
  array_name##_resize(array, 0); \
} \
\
static inline void array_name##_free(array_name##_t* array) \
{ \
  array_name##_clear(array); \
  if (array->dtors != NULL) \
    polymec_free(array->dtors); \
  if (array->data != NULL) \
    polymec_free(array->data); \
  polymec_free(array); \
} \
\
static inline bool array_name##_empty(array_name##_t* array) \
{ \
  return (array->size == 0); \
} \
\
static inline element* array_name##_find(array_name##_t* array, element value, array_name##_comparator comparator) \
{ \
  for (size_t i = 0; i < array->size; ++i) \
  { \
    if (comparator(value, array->data[i]) == 0) \
      return &array->data[i]; \
  } \
  return NULL; \
} \
static inline void array_name##_append_with_dtor(array_name##_t* array, element value, array_name##_dtor dtor) \
{ \
  array_name##_resize(array, array->size + 1); \
  array->data[array->size-1] = value; \
  if (dtor != NULL) \
  { \
    if (array->dtors == NULL) \
    { \
      array->dtors = (array_name##_dtor*)polymec_malloc(sizeof(array_name##_dtor) * array->capacity); \
      memset(array->dtors, 0, sizeof(array_name##_dtor) * array->capacity); \
    } \
    else \
    { \
      array->dtors = (array_name##_dtor*)polymec_realloc(array->dtors, sizeof(array_name##_dtor) * array->capacity); \
      memset(&array->dtors[array->size], 0, sizeof(array_name##_dtor) * (array->capacity - array->size)); \
    } \
    array->dtors[array->size-1] = dtor; \
  } \
} \
\
static inline void array_name##_append(array_name##_t* array, element value) \
{ \
  array_name##_append_with_dtor(array, value, NULL); \
} \
\
static inline void array_name##_insert_with_dtor(array_name##_t* array, size_t i, element value, array_name##_dtor dtor) \
{ \
  ASSERT(i <= array->size); \
  array_name##_resize(array, array->size + 1); \
  if (i < array->size-1) \
    memmove(&array->data[i+1], &array->data[i], sizeof(element) * (array->size-i)); \
  array->data[i] = value; \
  if (dtor != NULL) \
  { \
    if (array->dtors == NULL) \
    { \
      array->dtors = (array_name##_dtor*)polymec_malloc(sizeof(array_name##_dtor) * array->capacity); \
      memset(array->dtors, 0, sizeof(array_name##_dtor) * array->capacity); \
    } \
    else \
    { \
      array->dtors = (array_name##_dtor*)polymec_realloc(array->dtors, sizeof(array_name##_dtor) * array->capacity); \
      memmove(&array->dtors[i+1], &array->dtors[i], sizeof(array_name##_dtor) * (array->size-i)); \
    } \
    array->dtors[i] = dtor; \
  } \
} \
\
static inline void array_name##_insert(array_name##_t* array, size_t i, element value) \
{ \
  array_name##_insert_with_dtor(array, i, value, NULL); \
} \
\
static inline void array_name##_assign_with_dtor(array_name##_t* array, size_t i, element value, array_name##_dtor dtor) \
{ \
  ASSERT(i < array->size); \
  if (array->dtors[i] != NULL) \
    array->dtors[i](array->data[i]); \
  array->data[i] = value; \
  array->dtors[i] = dtor; \
} \
\
static inline void array_name##_assign(array_name##_t* array, size_t i, element value) \
{ \
  array_name##_assign_with_dtor(array, i, value, NULL); \
} \
\
static inline bool array_name##_next(array_name##_t* array, int* pos, element* value) \
{ \
  if (*pos < array->size) \
  { \
    *value = array->data[*pos]; \
    ++(*pos); \
    return true; \
  } \
  else \
    return false; \
} \
\
static inline void array_name##_release_data_and_free(array_name##_t* array) \
{ \
  array->data = NULL; \
  array_name##_free(array); \
} \

// Define some basic array types.
DEFINE_ARRAY(byte_array, uint8_t)
DEFINE_ARRAY(int_array, int)
DEFINE_ARRAY(index_array, index_t)
DEFINE_ARRAY(real_array, real_t)
DEFINE_ARRAY(string_array, char*)
DEFINE_ARRAY(ptr_array, void*)

#endif
