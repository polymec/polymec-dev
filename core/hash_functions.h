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

#ifndef POLYMEC_HASH_FUNCTIONS_H
#define POLYMEC_HASH_FUNCTIONS_H

#include <string.h>

// A solid string hash function by Dan Bernstein.
static inline int djb2_hash(unsigned char* str, int len)
{
  int hash = 5381;

  for (int i = 0; i < len; ++i)
  {
    int c = str[i];
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  }

  return hash;
}

// A slightly improved version of djb2 using xor.
static inline int djb2_xor_hash(unsigned char* str, int len)
{
  int hash = 5381;

  for (int i = 0; i < len; ++i)
  {
    int c = str[i];
    hash = ((hash << 5) + hash) ^ c;
  }

  return hash;
}

static inline int int_hash(int i)
{
  return djb2_xor_hash((unsigned char*)&i, sizeof(int));
}

static inline int int_pair_hash(int* i)
{
  return djb2_xor_hash((unsigned char*)i, 2*sizeof(int));
}

static inline int string_hash(char* str)
{
  return djb2_xor_hash((unsigned char*)str, strlen(str));
}

static inline int pointer_hash(void* p)
{
  return djb2_xor_hash((unsigned char*)p, sizeof(void*));
}

#endif
