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

#ifndef POLYMEC_COMPARATORS_H
#define POLYMEC_COMPARATORS_H

// These comparators return -1 if x < y, 0 if x == y, and 1 if x > y.

static inline int int_cmp(int x, int y)
{
  return (x < y) ? -1 : (x == y) ? 0 : 1;
}

static inline int double_cmp(double x, double y)
{
  return (x < y) ? -1 : (x == y) ? 0 : 1;
}

static inline bool int_equals(int x, int y)
{
  return (x == y);
}

static inline bool int_pair_equals(int* x, int* y)
{
  return ((x[0] == y[0]) && (x[1] == y[1]));
}

static inline bool string_equals(char* x, char* y)
{
  return (strcmp(x, y) == 0);
}

static inline bool pointer_equals(void* x, void* y)
{
  return (x == y);
}

// Not the right place for this guy, but it doesn't matter.
static inline void string_free(char* str)
{
  free(str);
}

#endif
