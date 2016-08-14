// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_COMPARATORS_H
#define POLYMEC_COMPARATORS_H

#include "core/polymec.h"

// These comparators return -1 if x < y, 0 if x == y, and 1 if x > y.

static inline int int_cmp(int x, int y)
{
  return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static inline int index_cmp(index_t x, index_t y)
{
  return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static inline int double_cmp(double x, double y)
{
  return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static inline int float_cmp(float x, float y)
{
  return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static inline int real_cmp(real_t x, real_t y)
{
  return (x < y) ? -1 : (x > y) ? 1 : 0;
}

static inline bool int_equals(int x, int y)
{
  return (x == y);
}

static inline bool index_equals(index_t x, index_t y)
{
  return (x == y);
}

static inline bool int_pair_equals(int* x, int* y)
{
  return ((x[0] == y[0]) && (x[1] == y[1]));
}

static inline bool index_pair_equals(index_t* x, index_t* y)
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

#endif
