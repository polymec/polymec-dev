#ifndef ARBI_VECTOR_H
#define ARBI_VECTOR_H

#include "core/arbi.h"
#include "arena/proto.h"

// A vector is a managed array that can be dynamically resized at your 
// convenience. Calling DEFINE_VECTOR_T(x) for a type x produces code for 
// a type called x_vector_t, and defines the following data structure and 
// interface:
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
// x_vector_foreach(vec, v) - Executes v(e) for each element e in vec.
//
// Other curiosities:
// x_vector_new_with_arena(arena, N) - Creates a new x_vector_t with N 
//                                     elements using the given arena.

#define DEFINE_VECTOR_T(element) \
typedef struct \
{ \
  element* data; \
  int size, capacity; \
  ARENA* arena; \
  bool close_arena; \
} element##_vector_t; \
\
typedef void (*element##_vector_visitor_t)(element e); \
static inline element##_vector_t* element##_vector_new_with_arena(ARENA* arena, int size) \
{ \
  ASSERT(size >= 0); \
  element##_vector_t* v = arena_malloc(arena, sizeof(element##_vector_t), 0); \
  v->arena = arena; \
  v->size = size; \
  v->capacity = 1; \
  while (v->capacity < v->size) \
    v->capacity *= 2; \
  v->data = arena_malloc(arena, sizeof(element)*v->capacity, 0); \
  v->close_arena = false; \
  return v; \
} \
\
static inline element##_vector_t* element##_vector_new(int size) \
{ \
  ARENA* a = arena_open(&arena_defaults, 0); \
  element##_vector_t* v = element##_vector_new_with_arena(a, size); \
  v->close_arena = true; \
  return v; \
} \
\
static inline void element##_vector_free(element##_vector_t* v) \
{ \
  arena_free(v->arena, v->data); \
  ARENA* arena = v->arena; \
  bool close_arena = v->close_arena; \
  arena_free(arena, v); \
  if (close_arena) \
    arena_close(arena); \
} \
\
static inline void element##_vector_reserve(element##_vector_t* v, int capacity) \
{ \
  int old_capacity = v->capacity; \
  while (capacity > v->capacity) \
    v->capacity *= 2; \
  if (v->capacity > old_capacity) \
    v->data = arena_realloc(v->arena, v->data, sizeof(element)*v->capacity, 0); \
} \
\
static inline void element##_vector_resize(element##_vector_t* v, int size) \
{ \
  if (size > v->size) \
    element##_vector_reserve(v, size); \
  v->size = size; \
} \
\
static inline void element##_vector_append(element##_vector_t* v, element e) \
{ \
  if (v->size == v->capacity) \
    element##_vector_reserve(v, v->size+1); \
  v->data[v->size] = e; \
  v->size++; \
} \
static inline void element##_vector_foreach(element##_vector_t* v, element##_vector_visitor_t visitor) \
{ \
  for (int i = 0; i < v->size; ++i) \
    visitor(v->data[i]); \
} \

#ifdef __cplusplus
extern "C" {
#endif

// Define some vectors.
DEFINE_VECTOR_T(int)
DEFINE_VECTOR_T(double)

#ifdef __cplusplus
}
#endif

#endif
