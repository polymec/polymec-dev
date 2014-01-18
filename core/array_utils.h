// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#ifndef POLYMEC_ARRAY_UTILS_H
#define POLYMEC_ARRAY_UTILS_H

#include "core/polymec.h"

// Integer binary search comparison function.
static inline int int_bsearch_comp(const void* l, const void* r)
{
  int li = *((int*)l), ri = *((int*)r);
  return (li < ri) ? -1
                   : (li > ri) ? 1
                               : 0;
}

// Real-valued binary search comparison function.
static inline int real_bsearch_comp(const void* l, const void* r)
{
  real_t li = *((real_t*)l), ri = *((real_t*)r);
  return (li < ri) ? -1
                   : (li > ri) ? 1
                               : 0;
}

// Executes a binary search for an element in a sorted array of integers, 
// returning a pointer to the element if it's found and NULL if it's not.
int* int_bsearch(int* array, int length, int element);

// Returns the index of the (sorted) array at which the desired element 
// appears (if it is present), or would appear (if it is not).
int int_lower_bound(int* array, int length, int element);

// Sorts (in-place) the elements in an array of integers. Uses qsort().
void int_qsort(int* array, int length);

// Executes a binary search for an element in a sorted array of real numbers, 
// returning a pointer to the element if it's found and NULL if it's not.
real_t* real_bsearch(real_t* array, int length, int element);

// Returns the index of the (sorted) array at which the desired element 
// appears (if it is present), or would appear (if it is not).
int real_lower_bound(real_t* array, int length, int element);

// Sorts (in-place) the elements in an array of reals. Uses qsort().
void real_qsort(real_t* array, int length);

#endif
