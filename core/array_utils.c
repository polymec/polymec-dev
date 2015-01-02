// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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

