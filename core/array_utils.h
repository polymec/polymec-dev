// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ARRAY_UTILS_H
#define POLYMEC_ARRAY_UTILS_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// Integer binary search comparison function.
static inline int int_bsearch_comp(const void* l, const void* r)
{
  int li = *((int*)l), ri = *((int*)r);
  return (li < ri) ? -1
                   : (li > ri) ? 1
                               : 0;
}

/// Index (large integer) binary search comparison function.
static inline int index_bsearch_comp(const void* l, const void* r)
{
  index_t li = *((index_t*)l), ri = *((index_t*)r);
  return (li < ri) ? -1
                   : (li > ri) ? 1
                               : 0;
}

/// Real-valued binary search comparison function.
static inline int real_bsearch_comp(const void* l, const void* r)
{
  real_t li = *((real_t*)l), ri = *((real_t*)r);
  return (li < ri) ? -1
                   : (li > ri) ? 1
                               : 0;
}

/// Integer pair binary search comparison function.
static inline int int_pair_bsearch_comp(const void* l, const void* r)
{
  int *li = (int*)l, *ri = (int*)r;
  return (li[0] < ri[0]) ? -1
                   : (li[0] == ri[0]) ? ((li[1] < ri[1]) ? -1 : (li[1] > ri[1]) ? 1 : 0)
                   : 1;
}

/// String-valued binary search comparison function.
static inline int string_bsearch_comp(const void* l, const void* r)
{
  char* li = *((char**)l), *ri = *((char**)r);
  return strcmp((const char*)li, (const char*)ri);
}

/// Fills an array of integers with the given value.
void int_fill(int* array, size_t length, int value);

/// Executes a linear search for an element in an unsorted array of integers,
/// returning a pointer to the element if it's found and NULL if it's not.
int* int_lsearch(int* array, size_t length, int element);

/// Executes a binary search for an element in a sorted array of integers,
/// returning a pointer to the element if it's found and NULL if it's not.
int* int_bsearch(int* array, size_t length, int element);

/// Returns the index of the (sorted) array at which the desired element
/// appears (if it is present), or would appear (if it is not).
size_t int_lower_bound(int* array, size_t length, int element);

/// Sorts (in-place) the elements in an array of integers. Uses qsort().
void int_qsort(int* array, size_t length);

/// Sorts the elements in an array of integers into the permutation array perm.
void int_qsort_to_perm(int* array, size_t length, size_t* perm);

/// Sorts (in-place) the elements in an array of integer pairs. Uses qsort().
void int_pair_qsort(int* array, size_t length);

/// Fills an array of indices with the given value.
void index_fill(index_t* array, size_t length, index_t value);

/// Executes a linear search for an element in an unsorted array of indices,
/// returning a pointer to the element if it's found and NULL if it's not.
index_t* index_lsearch(index_t* array, size_t length, index_t element);

/// Executes a binary search for an element in a sorted array of indices,
/// returning a pointer to the element if it's found and NULL if it's not.
index_t* index_bsearch(index_t* array, size_t length, index_t element);

/// Returns the index of the (sorted) array at which the desired element
/// appears (if it is present), or would appear (if it is not).
size_t index_lower_bound(index_t* array, size_t length, index_t element);

/// Sorts (in-place) the elements in an array of indices. Uses qsort().
void index_qsort(index_t* array, size_t length);

/// Sorts the elements in an array of indices into the permutation array perm.
void index_qsort_to_perm(index_t* array, size_t length, size_t* perm);

/// Fills an array of real numbers with the given value.
void real_fill(real_t* array, size_t length, real_t value);

/// Executes a linear search for an element in an unsorted array of real numbers,
/// returning a pointer to the element if it's found and NULL if it's not.
real_t* real_lsearch(real_t* array, size_t length, real_t element);

/// Executes a binary search for an element in a sorted array of real numbers,
/// returning a pointer to the element if it's found and NULL if it's not.
real_t* real_bsearch(real_t* array, size_t length, real_t element);

/// Returns the index of the (sorted) array at which the desired element
/// appears (if it is present), or would appear (if it is not).
size_t real_lower_bound(real_t* array, size_t length, real_t element);

/// Sorts (in-place) the elements in an array of reals. Uses qsort().
void real_qsort(real_t* array, size_t length);

/// Sorts the elements in an array of reals into the permutation array perm.
void real_qsort_to_perm(real_t* array, size_t length, size_t* perm);

///@}

#endif
