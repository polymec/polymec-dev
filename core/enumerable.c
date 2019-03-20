// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/enumerable.h"

bool enumerable_all(bool_array_t* values)
{
  bool all = true;
  for (size_t i = 0; i < values->size; ++i)
  {
    all = all && values->data[i];
    if (!all) break;
  }
  bool_array_free(values);
  return all;
}

bool enumerable_any(bool_array_t* values)
{
  bool any = false;
  for (size_t i = 0; i < values->size; ++i)
  {
    any = values->data[i];
    if (any) break;
  }
  bool_array_free(values);
  return any;
}

bool enumerable_none(bool_array_t* values)
{
  return !enumerable_any(values);
}

bool_array_t* enumerable_or(bool_array_t* values1, bool_array_t* values2)
{
  ASSERT(values1->size == values2->size);
  bool_array_t* result = bool_array_new_with_size(values1->size);
  for (size_t i = 0; i < result->size; ++i)
    result->data[i] = values1->data[i] || values2->data[i];
  bool_array_free(values1);
  bool_array_free(values2);
  return result;
}

bool_array_t* enumerable_and(bool_array_t* values1, bool_array_t* values2)
{
  ASSERT(values1->size == values2->size);
  bool_array_t* result = bool_array_new_with_size(values1->size);
  for (size_t i = 0; i < result->size; ++i)
    result->data[i] = values1->data[i] && values2->data[i];
  bool_array_free(values1);
  bool_array_free(values2);
  return result;
}

bool_array_t* enumerable_xor(bool_array_t* values1, bool_array_t* values2)
{
  ASSERT(values1->size == values2->size);
  bool_array_t* result = bool_array_new_with_size(values1->size);
  for (size_t i = 0; i < result->size; ++i)
  {
    result->data[i] = ((values1->data[i] && !values2->data[i]) ||
                       (!values1->data[i] && values2->data[i]));
  }
  bool_array_free(values1);
  bool_array_free(values2);
  return result;
}

bool_array_t* enumerable_not(bool_array_t* values)
{
  bool_array_t* result = bool_array_new_with_size(values->size);
  for (size_t i = 0; i < result->size; ++i)
    result->data[i] = !values->data[i];
  bool_array_free(values);
  return result;
}

// Generic function for defining compare_values implementations.
#define DEFINE_COMPARE_VALUES(datatype_prefix, datatype) \
bool_array_t* compare_##datatype_prefix##_values(datatype_prefix##_enumerable_generator_t* g1, \
                                                 datatype_prefix##_enumerable_generator_t* g2, \
                                                 bool (*compare)(datatype x, datatype y)) \
{ \
  bool_array_t* result = NULL; \
  if ((g1->array != NULL) && (g2->array != NULL)) \
  { \
    ASSERT(g1->num_values == g2->num_values); \
    result = bool_array_new_with_size(g1->num_values); \
    for (size_t i = 0; i < g1->num_values; ++i) \
    { \
      datatype x = g1->array[i]; \
      datatype y = g2->array[i]; \
      result->data[i] = compare(x, y); \
    } \
  } \
  else if ((g1->array != NULL) && (g2->array == NULL)) \
  { \
    result = bool_array_new_with_size(g1->num_values); \
    for (size_t i = 0; i < g1->num_values; ++i) \
    { \
      datatype x = g1->array[i]; \
      datatype y; \
      bool g2_has_val = g2->generate(g2->context, i, &y); \
      ASSERT(g2_has_val); \
      if (g2_has_val) \
        result->data[i] = compare(x, y); \
    } \
  } \
  else if ((g1->array == NULL) && (g2->array != NULL)) \
  { \
    result = bool_array_new_with_size(g2->num_values); \
    for (size_t i = 0; i < g2->num_values; ++i) \
    { \
      datatype x; \
      bool g1_has_val = g1->generate(g1->context, i, &x); \
      ASSERT(g1_has_val); \
      if (g1_has_val) \
      { \
        datatype y = g2->array[i]; \
        result->data[i] = compare(x, y); \
      } \
    } \
  } \
  else \
  { \
    result = bool_array_new(); \
    size_t i = 0; \
    while (true) \
    { \
      datatype x, y; \
      bool g1_has_val = g1->generate(g1->context, i, &x); \
      bool g2_has_val = g2->generate(g2->context, i, &y); \
      ASSERT((g1_has_val && g2_has_val) || \
             (!g1_has_val && !g2_has_val)); \
      if (!g1_has_val || !g2_has_val) break; \
      bool_array_append(result, compare(x, y)); \
      ++i; \
    } \
  } \
  release_ref(g1); \
  release_ref(g2); \
  return result; \
} \

DEFINE_COMPARE_VALUES(real, real_t)
DEFINE_COMPARE_VALUES(int, int)

// Generic function for defining functions for ordering operators.
#define DEFINE_ORDERING_CMP(cmp_name, datatype_prefix, datatype, cmp_func) \
bool_array_t* datatype_prefix##_##cmp_name(datatype_prefix##_enumerable_generator_t* g, \
                                           datatype value) \
{ \
  bool_array_t* result = NULL; \
  if (g->array != NULL) \
  { \
    result = bool_array_new_with_size(g->num_values); \
    for (size_t i = 0; i < g->num_values; ++i) \
    { \
      datatype x = g->array[i]; \
      result->data[i] = cmp_func(x, value); \
    } \
  } \
  else \
  { \
    result = bool_array_new(); \
    size_t i = 0; \
    while (true) \
    { \
      datatype x; \
      bool g_has_val = g->generate(g->context, i, &x); \
      if (!g_has_val) break; \
      bool_array_append(result, cmp_func(x, value)); \
      ++i; \
    } \
  } \
  release_ref(g); \
  return result; \
} \
\

static inline bool int_gt(int x, int y)
{
  return (x > y);
}

static inline bool int_gte(int x, int y)
{
  return (x >= y);
}

static inline bool int_lt(int x, int y)
{
  return (x < y);
}

static inline bool int_lte(int x, int y)
{
  return (x <= y);
}

static inline bool int_eq(int x, int y)
{
  return (x == y);
}

DEFINE_ORDERING_CMP(greater_than, int, int, int_gt)
DEFINE_ORDERING_CMP(greater_than_or_equal_to, int, int, int_gte)
DEFINE_ORDERING_CMP(less_than, int, int, int_lt)
DEFINE_ORDERING_CMP(less_than_or_equal_to, int, int, int_lte)
DEFINE_ORDERING_CMP(equal_to, int, int, int_eq)

static inline bool real_gt(real_t x, real_t y)
{
  return (x > y);
}

static inline bool real_gte(real_t x, real_t y)
{
  return (x >= y);
}

static inline bool real_lt(real_t x, real_t y)
{
  return (x < y);
}

static inline bool real_lte(real_t x, real_t y)
{
  return (x <= y);
}

DEFINE_ORDERING_CMP(greater_than, real, real_t, real_gt)
DEFINE_ORDERING_CMP(greater_than_or_equal_to, real, real_t, real_gte)
DEFINE_ORDERING_CMP(less_than, real, real_t, real_lt)
DEFINE_ORDERING_CMP(less_than_or_equal_to, real, real_t, real_lte)
DEFINE_ORDERING_CMP(equal_to, real, real_t, reals_equal)

// Generic function for printing values from generators
#define DEFINE_APPLY_TO_VALUES(datatype_prefix, datatype) \
void apply_to_##datatype_prefix##_values(void (*func)(void* context, datatype value), \
                                         void* context, \
                                         datatype_prefix##_enumerable_generator_t* g) \
{ \
  if (g->array != NULL) \
  { \
    for (size_t i = 0; i < g->num_values; ++i) \
    { \
      datatype x = g->array[i]; \
      func(context, x); \
    } \
  } \
  else \
  { \
    size_t i = 0; \
    while (true) \
    { \
      datatype x; \
      bool g_has_val = g->generate(g->context, i, &x); \
      if (!g_has_val) break; \
      func(context, x); \
      ++i; \
    } \
  } \
  release_ref(g); \
} \
\

DEFINE_APPLY_TO_VALUES(int, int)
DEFINE_APPLY_TO_VALUES(real, real_t)

