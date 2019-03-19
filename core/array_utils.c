// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include "core/array_utils.h"

// Sort to a permutation array.
typedef struct
{
  void* elem;
  size_t index;
  int (*cmp)(const void* l, const void* r);
} qsort_perm_t;

static int perm_cmp_func(const void* l, const void* r)
{
  const qsort_perm_t* ll = l;
  const qsort_perm_t* rr = r;
  return ll->cmp(ll->elem, rr->elem);
}

static void qsort_to_perm(void* array,
                          size_t length,
                          size_t elem_size,
                          int (*cmp)(const void* l, const void* r),
                          size_t* perm)
{
  qsort_perm_t elems[length];
  for (size_t i = 0; i < length; ++i)
  {
    elems[i].elem = &(((char*)array)[elem_size*i]);
    elems[i].index = i;
    elems[i].cmp = cmp;
  }
  qsort(elems, (size_t)length, sizeof(qsort_perm_t), perm_cmp_func);
  for (size_t i = 0; i < length; ++i)
    perm[i] = elems[i].index;
}

// This is the generic implementation of lower_bound().
static size_t lower_bound(void* array, size_t length, void* element, size_t elem_size, int (*comp)(const void*, const void*))
{
  size_t first = 0;
  size_t count = length, step;
  while (count > 0)
  {
    step = count/2;
    char* bytes = (char*)array;
    if (comp((void*)&bytes[elem_size*(first+step)], element) < 0)
    {
      first += step+1;
      count -= (step+1);
    }
    else
      count = step;
  }
  return first;
}

void int_fill(int* array, size_t length, int value)
{
  for (int i = 0; i < length; ++i)
    array[i] = value;
}

int* int_lsearch(int* array, size_t length, int element)
{
  for (int i = 0; i < length; ++i)
  {
    if (array[i] == element)
      return &array[i];
  }
  return NULL;
}

int* int_bsearch(int* array, size_t length, int element)
{
  return bsearch(&element, array, (size_t)length, sizeof(int), int_bsearch_comp);
}

size_t int_lower_bound(int* array, size_t length, int element)
{
  return lower_bound(array, length, &element, sizeof(int), int_bsearch_comp);
}

void int_qsort(int* array, size_t length)
{
  qsort(array, (size_t)length, sizeof(int), int_bsearch_comp);
}

void int_qsort_to_perm(int* array, size_t length, size_t* perm)
{
  qsort_to_perm(array, length, sizeof(int), int_bsearch_comp, perm);
}

void int_pair_qsort(int* array, size_t length)
{
  qsort(array, (size_t)length, 2*sizeof(int), int_pair_bsearch_comp);
}

void index_fill(index_t* array, size_t length, index_t value)
{
  for (int i = 0; i < length; ++i)
    array[i] = value;
}

index_t* index_lsearch(index_t* array, size_t length, index_t element)
{
  for (int i = 0; i < length; ++i)
  {
    if (array[i] == element)
      return &array[i];
  }
  return NULL;
}

index_t* index_bsearch(index_t* array, size_t length, index_t element)
{
  return bsearch(&element, array, (size_t)length, sizeof(index_t), index_bsearch_comp);
}

size_t index_lower_bound(index_t* array, size_t length, index_t element)
{
  return lower_bound(array, length, &element, sizeof(index_t), index_bsearch_comp);
}

void index_qsort(index_t* array, size_t length)
{
  qsort(array, (size_t)length, sizeof(index_t), index_bsearch_comp);
}

void index_qsort_to_perm(index_t* array, size_t length, size_t* perm)
{
  qsort_to_perm(array, length, sizeof(index_t), index_bsearch_comp, perm);
}

void real_fill(real_t* array, size_t length, real_t value)
{
  for (int i = 0; i < length; ++i)
    array[i] = value;
}

real_t* real_lsearch(real_t* array, size_t length, real_t element)
{
  for (int i = 0; i < length; ++i)
  {
    if (reals_equal(array[i], element))
      return &array[i];
  }
  return NULL;
}

real_t* real_bsearch(real_t* array, size_t length, real_t element)
{
  return bsearch(&element, array, (size_t)length, sizeof(real_t), real_bsearch_comp);
}

size_t real_lower_bound(real_t* array, size_t length, real_t element)
{
  return lower_bound(array, length, &element, sizeof(real_t), real_bsearch_comp);
}

void real_qsort(real_t* array, size_t length)
{
  qsort(array, (size_t)length, sizeof(real_t), real_bsearch_comp);
}

void real_qsort_to_perm(real_t* array, size_t length, size_t* perm)
{
  qsort_to_perm(array, length, sizeof(real_t), real_bsearch_comp, perm);
}

