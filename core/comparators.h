// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#ifndef POLYMEC_COMPARATORS_H
#define POLYMEC_COMPARATORS_H

// These comparators return -1 if x < y, 0 if x == y, and 1 if x > y.

static inline int int_cmp(int x, int y)
{
  return (x < y) ? -1 : (x == y) ? 0 : 1;
}

static inline int real_cmp(real_t x, real_t y)
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
