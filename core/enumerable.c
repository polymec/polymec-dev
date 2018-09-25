// Copyright (c) 2012-2018, Jeffrey N. Johnson
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

bool_array_t* compare_real_values(real_enumerable_generator_t* g1,
                                  real_enumerable_generator_t* g2,
                                  bool (*compare)(real_t x, real_t y))
{
  bool_array_t* result = NULL;
  if ((g1->array != NULL) && (g2->array != NULL))
  {
    ASSERT(g1->num_values == g2->num_values);
    result = bool_array_new_with_size(g1->num_values);
    for (size_t i = 0; i < g1->num_values; ++i)
    {
      real_t x = g1->array[i];
      real_t y = g2->array[i];
      result->data[i] = compare(x, y);
    }
  }
  else if ((g1->array != NULL) && (g2->array == NULL))
  {
    result = bool_array_new_with_size(g1->num_values);
    for (size_t i = 0; i < g1->num_values; ++i)
    {
      real_t x = g1->array[i];
      real_t y;
      bool g2_has_val = g2->generate(g2->context, i, &y);
      ASSERT(g2_has_val);
      result->data[i] = compare(x, y);
    }
  }
  else if ((g1->array == NULL) && (g2->array != NULL))
  {
    result = bool_array_new_with_size(g2->num_values);
    for (size_t i = 0; i < g2->num_values; ++i)
    {
      real_t x;
      bool g1_has_val = g1->generate(g1->context, i, &x);
      ASSERT(g1_has_val);
      real_t y = g2->array[i];
      result->data[i] = compare(x, y);
    }
  }
  else // if ((g1->array == NULL) && (g2->array == NULL))
  {
    result = bool_array_new();
    size_t i = 0;
    while (true)
    {
      real_t x, y;
      bool g1_has_val = g1->generate(g1->context, i, &x);
      bool g2_has_val = g2->generate(g2->context, i, &y);
      ASSERT((g1_has_val && g2_has_val) || 
             (!g1_has_val && !g2_has_val));
      if (!g1_has_val) break;
      bool_array_append(result, compare(x, y));
      ++i;
    }
  }

  // Let go of the generators.
  polymec_release(g1);
  polymec_release(g2);

  return result;
}

