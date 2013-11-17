// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_ARRAY_UTILS_H
#define POLYMEC_ARRAY_UTILS_H

// Column index comparison function for bsearch within supermatrix NRformat.
static inline int int_bsearch_comp(const void* l, const void* r)
{
  int li = *((int*)l), ri = *((int*)r);
  return (li < ri) ? -1
                   : (li > ri) ? 1
                               : 0;
}

// Executes a binary search for an element in a sorted array of integers, 
// returning a pointer to the element if it's found and NULL if it's not.
int* int_bsearch(int* array, int length, int element);

// Sorts (in-place) the elements in an array of integers. Uses qsort().
void int_qsort(int* array, int length);

#endif
