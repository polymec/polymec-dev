// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ENUMERABLE_H
#define POLYMEC_ENUMERABLE_H

#include "core/array.h"

/// \addtogroup core core
///@{

/// Returns true if all of the given values are true, false otherwise.
/// The array of values is consumed.
bool enumerable_all(bool_array_t* values);

/// Returns true if any of the given values are true, false if none are.
/// The array of values is consumed.
bool enumerable_any(bool_array_t* values);

/// Returns true if none of the given values are true, false if any are.
/// The array of values is consumed.
bool enumerable_none(bool_array_t* values);

/// \macro ALL
/// Evaluates to true if all of the given values in the expression are true, 
/// and to false otherwise.
#define ALL(x) enumerable_all(x)

/// \macro ANY
/// Evaluates to true if any of the given values in the expression are true, 
/// and to false if none are true.
#define ANY(x) enumerable_any(x)

/// \macro NONE
/// Evaluates to true if none of the given values in the expression are true, 
/// and to false if any are true.
#define NONE(x) enumerable_none(x)

/// \def DEFINE_ENUMERABLE_GENERATOR
/// Defines a generator class for various enumerable value types.
#define DEFINE_ENUMERABLE_GENERATOR(generator_name, value) \
typedef value generator_name##_value_t; \
typedef struct generator_name##_t generator_name##_t; \
struct generator_name##_t \
{ \
  size_t position; \
	value* array; \
	size_t num_values; \
  bool (*generate)(void* context, size_t index, value* val); \
  void* context; \
  void (*dtor)(void* context); \
}; \
\
static inline void generator_name##_free(void* ctx) \
{ \
  generator_name##_t* gen = (generator_name##_t*)ctx; \
  if (gen->dtor != NULL) \
  { \
    if (gen->array != NULL) \
      gen->dtor((void*)(gen->array)); \
    else \
      gen->dtor(gen->context); \
  } \
} \
\
static inline generator_name##_t* generator_name##_from_array(value* values, size_t num_values, bool assume_ownership) \
{ \
  generator_name##_t* gen = (generator_name##_t*)polymec_gc_malloc(sizeof(generator_name##_t), generator_name##_free); \
  gen->position = 0; \
  gen->array = values; \
  gen->num_values = num_values; \
  gen->generate = NULL; \
  gen->context = NULL; \
  if (assume_ownership) \
    gen->dtor = polymec_free; \
  else \
    gen->dtor = NULL; \
  return gen; \
} \
\
static inline generator_name##_t* generator_name##_new(bool (*generate)(void* context, size_t index, value* val), \
                                                       void* context, \
                                                       void (*dtor)(void*)) \
{ \
  generator_name##_t* gen = (generator_name##_t*)polymec_malloc(sizeof(generator_name##_t)); \
  gen->position = 0; \
  gen->array = NULL; \
  gen->num_values = 0; \
  gen->generate = generate; \
  gen->context = context; \
  gen->dtor = dtor; \
  return gen; \
} \

// Define some enumerable generators.
DEFINE_ENUMERABLE_GENERATOR(real_enumerable_generator, real_t)

// Fancy type-generic macros aren't available to C++.
#ifndef __cplusplus

/// Compare the values produced by two enumerable generators x and y.
#define compare_values(x, y, cmp) _Generic((x), \
                       real_enumerable_generator_t*: compare_real_values)(x, y, cmp)

#endif // ifndef __cplusplus

bool_array_t* compare_real_values(real_enumerable_generator_t* g1,
                                  real_enumerable_generator_t* g2,
                                  bool (*compare)(real_t x, real_t y));


///@}


#endif
