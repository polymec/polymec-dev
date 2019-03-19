// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

/// \addtogroup enumerable enumerable
/// Enumerable objects are objects that contain values that can be ordered
/// and enumerated. That is: the values in an enumerable object can be mapped
/// to the positive integers. This allows us to analyze and compare values
/// in enumerable objects in high-level ways.
///
/// For example, we can use expressions that represent
/// * whether any value in an enumeration is negative:
///       \ref ANY(\ref less_than(enumerable, 0.0))
/// * whether all (real) values in two enumerable objects are equal:
///       \ref ALL(\ref compare_values(first_enumerable, second_enumerable, \ref reals_equal))
///
/// and so on.
///
/// All of the `enumerable` arguments above are represented by objects called
/// _enumerable generators_. You can define a specific enumeration for any
/// data type by constructing an `enumerable_generator` from its data and
/// passing it to one of the enumerable-friendly functions in this file.
///@{

/// Returns true if all of the given values are true, false otherwise.
/// The array of values is consumed.
/// \memberof enumerable
bool enumerable_all(bool_array_t* values);

/// Returns true if any of the given values are true, false if none are.
/// The array of values is consumed.
/// \memberof enumerable
bool enumerable_any(bool_array_t* values);

/// Returns true if none of the given values are true, false if any are.
/// The array of values is consumed.
/// \memberof enumerable
bool enumerable_none(bool_array_t* values);

/// Returns a new boolean array whose values are the pairwise OR of the values
/// of the two argument arrays. Both argument arrays are consumed.
/// \memberof enumerable
bool_array_t* enumerable_or(bool_array_t* values1, bool_array_t* values2);

/// Returns a new boolean array whose values are the pairwise AND of the values
/// of the two argument arrays. Both argument arrays are consumed.
/// \memberof enumerable
bool_array_t* enumerable_and(bool_array_t* values1, bool_array_t* values2);

/// Returns a new boolean array whose values are the pairwise XOR of the values
/// of the two argument arrays. Both argument arrays are consumed.
/// \memberof enumerable
bool_array_t* enumerable_xor(bool_array_t* values1, bool_array_t* values2);

/// Returns a new boolean array whose values are the negation of the values
/// of the argument array. The argument array is consumed.
/// \memberof enumerable
bool_array_t* enumerable_not(bool_array_t* values);

/// \macro ALL
/// Evaluates to true if all of the given values in the expression are true,
/// and to false otherwise.
/// \memberof enumerable
#define ALL(x) enumerable_all(x)

/// \macro ANY
/// Evaluates to true if any of the given values in the expression are true,
/// and to false if none are true.
/// \memberof enumerable
#define ANY(x) enumerable_any(x)

/// \macro NONE
/// Evaluates to true if none of the given values in the expression are true,
/// and to false if any are true.
/// \memberof enumerable
#define NONE(x) enumerable_none(x)

/// \macro OR
/// Given two boolean arrays, produces a boolean array whose values are the
/// pairwise OR of the originals.
/// \memberof enumerable
#define OR(x, y) enumerable_or(x, y)

/// \macro AND
/// Given two boolean arrays, produces a boolean array whose values are the
/// pairwise AND of the originals.
#define AND(x, y) enumerable_and(x, y)

/// \macro XOR
/// Given two boolean arrays, produces a boolean array whose values are the
/// pairwise XOR of the originals.
#define XOR(x, y) enumerable_xor(x, y)

/// \macro NOT
/// Produces a boolean array whose values are the negation of those of the
/// argument.
#define NOT(x) enumerable_not(x)

/// \def DEFINE_ENUMERABLE_GENERATOR
/// Defines a generator class for various enumerable value types.
/// \refcounted
/// \memberof enumerable
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
  generator_name##_t* gen = (generator_name##_t*)polymec_refcounted_malloc(sizeof(generator_name##_t), generator_name##_free); \
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
DEFINE_ENUMERABLE_GENERATOR(int_enumerable_generator, int)
DEFINE_ENUMERABLE_GENERATOR(real_enumerable_generator, real_t)

// Fancy type-generic macros aren't available to C++.
#ifndef __cplusplus

/// Compares the values produced by two enumerable generators x and y, using
/// a pairwise comparator cmp(x, y) that returns true or false.
/// \memberof enumerable
#define compare_values(x, y, cmp) _Generic((x), \
                       int_enumerable_generator_t*: compare_int_values, \
                       real_enumerable_generator_t*: compare_real_values)(x, y, cmp)

/// Determines whether the values produced by an enumerable generator x are
/// greater than the given fixed value.
/// \memberof enumerable
#define greater_than(x, y, cmp) _Generic((x), \
                     int_enumerable_generator_t*: int_greater_than, \
                     real_enumerable_generator_t*: real_greater_than)(x, y)

/// Determines whether the values produced by an enumerable generator x are
/// greater than or equal to the given fixed value.
/// \memberof enumerable
#define greater_than_or_equal_to(x, y, cmp) _Generic((x), \
                     int_enumerable_generator_t*: int_greater_than_or_equal_to, \
                     real_enumerable_generator_t*: real_greater_than_or_equal_to)(x, y)

/// Determines whether the values produced by an enumerable generator x are
/// less than the given fixed value.
/// \memberof enumerable
#define less_than(x, y, cmp) _Generic((x), \
                     int_enumerable_generator_t*: int_less_than, \
                     real_enumerable_generator_t*: real_less_than)(x, y)

/// Determines whether the values produced by an enumerable generator x are
/// less than or equal to the given fixed value.
/// \memberof enumerable
#define less_than_or_equal_to(x, y, cmp) _Generic((x), \
                     int_enumerable_generator_t*: int_less_than_or_equal_to, \
                     real_enumerable_generator_t*: real_less_than_or_equal_to)(x, y)

/// Determines whether the values produced by an enumerable generator x are
/// equal to the given fixed value.
/// \memberof enumerable
#define equal_to(x, y, cmp) _Generic((x), \
                     int_enumerable_generator_t*: int_equal_to, \
                     real_enumerable_generator_t*: real_equal_to)(x, y)

/// Applies the given function (with the given supplied context pointer) to each value
/// produced by an enumerable generator.
/// \param func A function that takes a context pointer and a value
/// \memberof enumerable
#define apply_to_values(func, context, x) _Generic((x), \
                     int_enumerable_generator_t*: apply_to_int_values, \
                     real_enumerable_generator_t*: apply_to_real_values)(func, context, x)

#endif // ifndef __cplusplus

// compare_values
bool_array_t* compare_int_values(int_enumerable_generator_t* g1,
                                 int_enumerable_generator_t* g2,
                                 bool (*compare)(int x, int y));
bool_array_t* compare_real_values(real_enumerable_generator_t* g1,
                                  real_enumerable_generator_t* g2,
                                  bool (*compare)(real_t x, real_t y));

// Ordering comparison functions
bool_array_t* int_greater_than(int_enumerable_generator_t* g, int value);
bool_array_t* int_greater_than_or_equal_to(int_enumerable_generator_t* g, int value);
bool_array_t* int_less_than(int_enumerable_generator_t* g, int value);
bool_array_t* int_less_than_or_equal_to(int_enumerable_generator_t* g, int value);
bool_array_t* int_equal_to(int_enumerable_generator_t* g, int value);

bool_array_t* real_greater_than(real_enumerable_generator_t* g, real_t value);
bool_array_t* real_greater_than_or_equal_to(real_enumerable_generator_t* g, real_t value);
bool_array_t* real_less_than(real_enumerable_generator_t* g, real_t value);
bool_array_t* real_less_than_or_equal_to(real_enumerable_generator_t* g, real_t value);
bool_array_t* real_equal_to(real_enumerable_generator_t* g, real_t value);

// apply_to_values
void apply_to_int_values(void (*func)(void* context, int value),
                         void* context,
                         int_enumerable_generator_t* g);
void apply_to_real_values(void (*func)(void* context, real_t value),
                          void* context,
                          real_enumerable_generator_t* g);

//
///@}


#endif
