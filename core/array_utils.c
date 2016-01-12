// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include "core/array_utils.h"

// This is the generic implementation of lower_bound(). 
static int lower_bound(void* array, int length, void* element, size_t elem_size, int (*comp)(const void*, const void*))
{
  int first = 0;
  int count = length, step;
  while (count > 0)
  {
    step = count/2;
    char* bytes = (char*)array;
    if (comp((void*)&bytes[elem_size*(first+step)], element) < 0)
    {
      first += (step+1);
      count -= (step+1);
    }
    else
      count = step;
  }
  return first;
}

int* int_lsearch(int* array, int length, int element)
{
  for (int i = 0; i < length; ++i)
  {
    if (array[i] == element)
      return &array[i];
  }
  return NULL;
}

int* int_bsearch(int* array, int length, int element)
{
  return bsearch(&element, array, (size_t)length, sizeof(int), int_bsearch_comp);
}

int int_lower_bound(int* array, int length, int element)
{
  return lower_bound(array, length, &element, sizeof(int), int_bsearch_comp);
}

void int_qsort(int* array, int length)
{
  qsort(array, (size_t)length, sizeof(int), int_bsearch_comp);
}

void int_pair_qsort(int* array, int length)
{
  qsort(array, (size_t)length, 2*sizeof(int), int_pair_bsearch_comp);
}

index_t* index_lsearch(index_t* array, int length, index_t element)
{
  for (int i = 0; i < length; ++i)
  {
    if (array[i] == element)
      return &array[i];
  }
  return NULL;
}

index_t* index_bsearch(index_t* array, int length, index_t element)
{
  return bsearch(&element, array, (size_t)length, sizeof(index_t), index_bsearch_comp);
}

int index_lower_bound(index_t* array, int length, index_t element)
{
  return lower_bound(array, length, &element, sizeof(index_t), index_bsearch_comp);
}

void index_qsort(index_t* array, int length)
{
  qsort(array, (size_t)length, sizeof(index_t), index_bsearch_comp);
}

real_t* real_lsearch(real_t* array, int length, real_t element)
{
  for (int i = 0; i < length; ++i)
  {
    if (array[i] == element)
      return &array[i];
  }
  return NULL;
}

real_t* real_bsearch(real_t* array, int length, real_t element)
{
  return bsearch(&element, array, (size_t)length, sizeof(real_t), real_bsearch_comp);
}

int real_lower_bound(real_t* array, int length, real_t element)
{
  return lower_bound(array, length, &element, sizeof(real_t), real_bsearch_comp);
}

void real_qsort(real_t* array, int length)
{
  qsort(array, (size_t)length, sizeof(real_t), real_bsearch_comp);
}

