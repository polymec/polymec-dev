#ifndef POLYMEC_VECTOR_H
#define POLYMEC_VECTOR_H

#include <stdlib.h>
#include "core/polymec.h"
#include "arena/proto.h"

// A vector is a managed array that can be dynamically resized at your 
// convenience. Calling DEFINE_VECTOR(x_vector, x) for a type x produces 
// code for a type called x_vector_t, and defines the following data 
// structure and interface:
//
// Suppose we have an x_vector_t* called vec. Then the following fields 
// are defined:
//
// vec->data     - An array of x.
// vec->size     - The number of x elements in the array.
//                 (This field must not be modified manually.)
// vec->capacity - The maximum number of x that can presently be stored.
//                 (This field must not be modified manually.)
// vec->arena,
// vec->close_arena - Used internally for memory allocation.
//
// Basic interface:
//
// vec = x_vector_new(N)    - Creates a new x_vector_t with N elements.
// x_vector_free(vec)       - Destroys vec, freeing its memory.
// x_vector_reserve(vec, N) - Reserves space for N elements within vec.
//                            Does not alter vec->size.
// x_vector_resize(vec, N)  - Resizes vec space for N elements within vec.
// x_vector_append(vec, e)  - Appends e to the end of vec, resizing it.
//
// Other curiosities:
// x_vector_new_with_arena(arena, N) - Creates a new x_vector_t with N 
//                                     elements using the given arena.

#define DEFINE_VECTOR(vector_name, element) \
typedef struct \
{ \
  element* data; \
  int size, capacity; \
  ARENA* arena; \
  bool close_arena; \
} vector_name##_t; \
\
typedef void (*vector_name##_visitor)(element e); \
static inline vector_name##_t* vector_name##_new_with_arena(ARENA* arena, int size) \
{ \
  ASSERT(size >= 0); \
  vector_name##_t* v = arena_malloc(arena, sizeof(vector_name##_t), 0); \
  v->arena = arena; \
  v->size = size; \
  v->capacity = 1; \
  while (v->capacity < v->size) \
    v->capacity *= 2; \
  v->data = arena_malloc(arena, sizeof(element)*v->capacity, 0); \
  return v; \
} \
\
static inline vector_name##_t* vector_name##_new(int size) \
{ \
  ASSERT(size >= 0); \
  vector_name##_t* v = malloc(sizeof(vector_name##_t)); \
  v->arena = NULL; \
  v->size = size; \
  v->capacity = 1; \
  while (v->capacity < v->size) \
    v->capacity *= 2; \
  v->data = malloc(sizeof(element)*v->capacity); \
  return v; \
} \
\
static inline void vector_name##_free(vector_name##_t* v) \
{ \
  if (v->arena != NULL) \
  { \
    arena_free(v->arena, v->data); \
    ARENA* arena = v->arena; \
    arena_free(arena, v); \
  } \
  else \
  { \
    free(v->data); \
    free(v); \
  } \
} \
\
static inline void vector_name##_reserve(vector_name##_t* v, int capacity) \
{ \
  int old_capacity = v->capacity; \
  while (capacity > v->capacity) \
    v->capacity *= 2; \
  if (v->capacity > old_capacity) \
  { \
    if (v->arena != NULL) \
      v->data = arena_realloc(v->arena, v->data, sizeof(element)*v->capacity, 0); \
    else \
      v->data = realloc(v->data, sizeof(element)*v->capacity); \
  } \
} \
\
static inline void vector_name##_resize(vector_name##_t* v, int size) \
{ \
  if (size > v->size) \
    vector_name##_reserve(v, size); \
  v->size = size; \
} \
\
static inline void vector_name##_append(vector_name##_t* v, element e) \
{ \
  if (v->size == v->capacity) \
    vector_name##_reserve(v, v->size+1); \
  v->data[v->size] = e; \
  v->size++; \
} \
\

#ifdef __cplusplus
extern "C" {
#endif

// Define some vectors.
DEFINE_VECTOR(int_vector, int)
DEFINE_VECTOR(double_vector, double)

#ifdef __cplusplus
}
#endif

#endif
